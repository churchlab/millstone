import datetime
import os
import pickle
import subprocess
import re

from Bio import SeqIO
from celery import group
from celery import task
from django.conf import settings

from genome_finish.graph_contig_placement import graph_contig_placement
from genome_finish.millstone_de_novo_fns import add_paired_mates
from genome_finish.millstone_de_novo_fns import filter_low_qual_read_pairs
from genome_finish.millstone_de_novo_fns import filter_out_unpaired_reads
from genome_finish.millstone_de_novo_fns import get_altalign_reads
from genome_finish.millstone_de_novo_fns import get_piled_reads
from genome_finish.millstone_de_novo_fns import get_clipped_reads_smart
from genome_finish.millstone_de_novo_fns import get_avg_genome_coverage
from genome_finish.millstone_de_novo_fns import get_unmapped_reads
from main.models import Contig
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.model_utils import get_dataset_with_type
from pipeline.read_alignment import get_insert_size_mean_and_stdev
from pipeline.read_alignment_util import ensure_bwa_index
from utils.bam_utils import concatenate_bams
from utils.bam_utils import index_bam
from utils.bam_utils import make_bam
from utils.bam_utils import make_sam
from utils.bam_utils import rmdup
from utils.bam_utils import sort_bam_by_coordinate
from utils.bam_utils import sort_bam_by_name
from utils.data_export_util import export_contig_list_as_vcf
from utils.data_export_util import export_var_dict_list_as_vcf
from utils.import_util import add_dataset_to_entity
from utils.import_util import prepare_ref_genome_related_datasets
from utils.jbrowse_util import add_bam_file_track
from variants.vcf_parser import parse_vcf

# Default args for velvet assembly
VELVETH_BINARY = settings.TOOLS_DIR + '/velvet/velveth'
VELVETG_BINARY = settings.TOOLS_DIR + '/velvet/velvetg'

VELVET_COVERAGE_CUTOFF = 30
VELVET_HASH_LENGTH = 21
VELVET_MIN_CONTIG_LENGTH = 200

# EXPERIMENTALLY BEST COVERAGE RATIOS
VELVET_CONTIG_COVERAGE_EXPECTED = 0.8
VELVET_CONTIG_COVERAGE_CUTOFF = 0.2

DEFAULT_VELVET_OPTS = {
    'velveth': {
        'hash_length': VELVET_HASH_LENGTH
    },
    'velvetg': {
        'read_trkg': 'yes',
        'cov_cutoff': VELVET_COVERAGE_CUTOFF,
        'min_contig_lgth': VELVET_MIN_CONTIG_LENGTH
    }
}

NUM_CONTIGS_TO_EVALUATE = 1000


def run_de_novo_assembly_pipeline(sample_alignment_list,
        sv_read_classes={}, input_velvet_opts={},
        overwrite=True):

    for sample_alignment in sample_alignment_list:
        sample_alignment.data['assembly_status'] = (
                ExperimentSampleToAlignment.ASSEMBLY_STATUS.QUEUED)
        sample_alignment.save()

    # Ensure reference genome fasta has bwa index
    ref_genome = sample_alignment_list[0].alignment_group.reference_genome
    ref_genome_fasta = ref_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    ensure_bwa_index(ref_genome_fasta)

    generate_contigs_async_task = generate_contigs_multi_sample(
            sample_alignment_list)
    async_result = generate_contigs_async_task.apply_async()

    return async_result


def kmer_coverage(C, L, k):
    """Converts contig coverage to kmer coverage

    Args:
        * C: Coverage
        * L: Template length
        * k: hash length
    """
    return C * (L - k + 1) / float(L)


def generate_contigs_multi_sample(sample_alignment_list):
    """Generate contigs for each ExperimentSampleToAlignment in
    the passed list in the order that they are displayed in the UI
    """
    generate_contigs_task_list = []
    for sample_alignment in sorted(sample_alignment_list,
            key=lambda x: x.experiment_sample.label):
        generate_contigs_task_list.append(
                generate_contigs.si(sample_alignment))

    task_signature = group(generate_contigs_task_list)
    for sample_alignment in sample_alignment_list:

        task_signature = task_signature | parse_variants_from_vcf.si(
                sample_alignment)

    return task_signature


@task
def generate_contigs(sample_alignment,
        sv_read_classes={}, input_velvet_opts={},
        overwrite=True):

    # Set assembly status for UI
    sample_alignment.data['assembly_status'] = (
                ExperimentSampleToAlignment.ASSEMBLY_STATUS.ASSEMBLING)
    sample_alignment.save()

    # Grab reference genome fasta path, ensure indexed
    reference_genome = sample_alignment.alignment_group.reference_genome
    ref_fasta_dataset = reference_genome.dataset_set.get_or_create(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA)[0]
    prepare_ref_genome_related_datasets(reference_genome, ref_fasta_dataset)

    # Make data_dir directory to house genome_finishing files
    assembly_dir = os.path.join(
            sample_alignment.get_model_data_dir(),
            'assembly')

    # Make assembly directory if it does not exist
    if not os.path.exists(assembly_dir):
        os.mkdir(assembly_dir)

    data_dir_counter = 0
    data_dir = os.path.join(assembly_dir, str(data_dir_counter))
    while(os.path.exists(data_dir)):
        data_dir_counter += 1
        data_dir = os.path.join(assembly_dir, str(data_dir_counter))
    os.mkdir(data_dir)

    # Get a bam of sorted SV indicants with pairs
    sv_indicants_bam = get_sv_indicating_reads(sample_alignment,
            sv_read_classes, overwrite=overwrite)

    prev_dataset = get_dataset_with_type(
            sample_alignment,
            Dataset.TYPE.BWA_FOR_DE_NOVO_ASSEMBLY)

    if overwrite and prev_dataset:
        prev_dataset.delete()

    if overwrite or prev_dataset is None:

        sv_indicants_sorted_bam = (os.path.splitext(sv_indicants_bam)[0] +
                '.coordinate_sorted.bam')

        # Bam needs to be coordinated sorted to index
        sort_bam_by_coordinate(sv_indicants_bam, sv_indicants_sorted_bam)

        # Bam needs to be indexed for jbrowse
        index_bam(sv_indicants_sorted_bam)

        for_assembly_dataset = add_dataset_to_entity(
                sample_alignment,
                Dataset.TYPE.BWA_FOR_DE_NOVO_ASSEMBLY,
                Dataset.TYPE.BWA_FOR_DE_NOVO_ASSEMBLY,
                filesystem_location=sv_indicants_sorted_bam)

        for_assembly_dataset.save()

        add_bam_file_track(reference_genome,
                sample_alignment, Dataset.TYPE.BWA_FOR_DE_NOVO_ASSEMBLY)

    velvet_opts = dict(DEFAULT_VELVET_OPTS)

    # Find insertion metrics
    ins_length, ins_length_sd = get_insert_size_mean_and_stdev(
            sample_alignment)
    velvet_opts['velvetg']['ins_length'] = ins_length
    velvet_opts['velvetg']['ins_length_sd'] = ins_length_sd

    # Find expected coverage
    avg_read_coverage = get_avg_genome_coverage(
            sample_alignment)

    # Calculate expected coverage in kmers
    genome_kmer_coverage = kmer_coverage(avg_read_coverage, ins_length,
            velvet_opts['velveth']['hash_length'])
    exp_cov = genome_kmer_coverage * VELVET_CONTIG_COVERAGE_EXPECTED
    velvet_opts['velvetg']['exp_cov'] = exp_cov

    # # Set cov cutoff
    cov_cutoff = genome_kmer_coverage * VELVET_CONTIG_COVERAGE_CUTOFF
    velvet_opts['velvetg']['cov_cutoff'] = cov_cutoff

    # Update velvet_opts with input_velvet_opts
    for shallow_key in ['velveth', 'velvetg']:
        if shallow_key in input_velvet_opts:
            for deep_key in input_velvet_opts[shallow_key]:
                velvet_opts[shallow_key][deep_key] = (
                        input_velvet_opts[shallow_key][deep_key])

    # Perform velvet assembly
    contig_files = assemble_with_velvet(
            data_dir, velvet_opts, sv_indicants_bam,
            sample_alignment)


    sample_alignment.data['assembly_status'] = (
            ExperimentSampleToAlignment.ASSEMBLY_STATUS.COMPLETED)
    sample_alignment.save()

    return contig_files


def get_sv_indicating_reads(sample_alignment, input_sv_indicant_classes={},
        overwrite=False):

    sv_indicant_keys = [
            Dataset.TYPE.BWA_ALTALIGN,
            Dataset.TYPE.BWA_PILED,
            Dataset.TYPE.BWA_CLIPPED,
            Dataset.TYPE.BWA_SPLIT,
            Dataset.TYPE.BWA_UNMAPPED,
            Dataset.TYPE.BWA_DISCORDANT
    ]

    sv_indicant_class_to_filename_suffix = {
            Dataset.TYPE.BWA_ALTALIGN: 'altalign',
            Dataset.TYPE.BWA_PILED: 'piled',
            Dataset.TYPE.BWA_CLIPPED: 'clipped',
            Dataset.TYPE.BWA_SPLIT: 'split',
            Dataset.TYPE.BWA_UNMAPPED: 'unmapped',
            Dataset.TYPE.BWA_DISCORDANT: 'discordant'
    }

    sv_indicant_class_to_generator = {
            Dataset.TYPE.BWA_ALTALIGN: get_altalign_reads,
            Dataset.TYPE.BWA_PILED: get_piled_reads,
            Dataset.TYPE.BWA_CLIPPED: lambda i, o: get_clipped_reads_smart(
                    i, o,
                    phred_encoding=sample_alignment.experiment_sample.data.get(
                            'phred_encoding', None)),
            Dataset.TYPE.BWA_UNMAPPED: lambda i, o: get_unmapped_reads(
                    i, o, avg_phred_cutoff=20)
    }

    default_sv_indicant_classes = {
            Dataset.TYPE.BWA_ALTALIGN: False,
            Dataset.TYPE.BWA_PILED: False,
            Dataset.TYPE.BWA_CLIPPED: True,
            Dataset.TYPE.BWA_SPLIT: True,
            Dataset.TYPE.BWA_UNMAPPED: True,
            Dataset.TYPE.BWA_DISCORDANT: True
    }
    default_sv_indicant_classes.update(input_sv_indicant_classes)

    # Grab alignment bam file-path
    alignment_bam = get_dataset_with_type(
            sample_alignment,
            Dataset.TYPE.BWA_ALIGN).get_absolute_location()

    # HACK: Filter out unpaired mates
    alignment_no_unpaired_bam = (os.path.splitext(alignment_bam)[0] +
            '.no_unpaired.bam')
    if not os.path.exists(alignment_no_unpaired_bam):
        filter_out_unpaired_reads(alignment_bam, alignment_no_unpaired_bam)

    # Use no unpaired filtered bam
    alignment_bam = alignment_no_unpaired_bam

    # Index it
    index_bam(alignment_bam)

    # Get SV indicating reads
    sv_bams_list = []
    alignment_file_prefix = os.path.join(
            sample_alignment.get_model_data_dir(),
            'bwa_align')

    # Helper function for getting sv read datasets
    def _get_or_create_sv_dataset(key):
        dataset_query = sample_alignment.dataset_set.filter(type=key)

        if dataset_query.exists() and not overwrite or (
                dataset_query.exists() and
                key not in sv_indicant_class_to_generator):
            assert len(dataset_query) == 1
            return dataset_query[0]
        elif dataset_query.exists() and overwrite and (
                key in sv_indicant_class_to_generator):
            assert len(dataset_query) == 1
            dataset_query[0].delete()

        if overwrite and key in sv_indicant_class_to_generator or (
                not dataset_query.exists()):
            dataset_path = '.'.join([
                    alignment_file_prefix,
                    sv_indicant_class_to_filename_suffix[key],
                    'bam'
                    ])
            generator = sv_indicant_class_to_generator[key]
            generator(alignment_bam, dataset_path)

            return add_dataset_to_entity(
                    sample_alignment,
                    key,
                    key,
                    filesystem_location=dataset_path)

    # Aggregate SV indicants
    for key in sv_indicant_keys:
        if default_sv_indicant_classes[key]:
            dataset = _get_or_create_sv_dataset(key)
            sv_bams_list.append(dataset.get_absolute_location())

    # Make some bam tracks for read classes
    jbrowse_classes = [Dataset.TYPE.BWA_DISCORDANT]
    reference_genome = sample_alignment.alignment_group.reference_genome
    for dataset_type in jbrowse_classes:
        bam_path = sample_alignment.dataset_set.get(
                type=dataset_type).get_absolute_location()
        index_bam(bam_path)
        add_bam_file_track(reference_genome,
        sample_alignment, dataset_type)

    # Create compilation filename prefix
    suffixes = [sv_indicant_class_to_filename_suffix[k]
            for k in sv_indicant_keys
            if default_sv_indicant_classes[k]]

    suffix_string = '_'.join(sorted(suffixes))
    compilation_prefix = '.'.join([
            alignment_file_prefix,
            suffix_string])

    SV_indicants_filtered = compilation_prefix + '.with_pairs.filtered.bam'
    if os.path.exists(SV_indicants_filtered) and not overwrite:
        print ('WARNING: Requested SV indicants bam file: ' +
                SV_indicants_filtered +
                ' already exists and will be returned by this function.  ' +
                'To overwrite this file pass the keyword overwrite=True')
        return SV_indicants_filtered
    if overwrite:
        print ('WARNING: overwrite is True, so SV read bam datasets ' +
                'are being overwritten')

    # Aggregate SV indicants
    print 'concatenating sv indicants'
    SV_indicants_bam = compilation_prefix + '.bam'
    concatenate_bams(
            sv_bams_list,
            SV_indicants_bam)

    # Remove duplicates
    print 'removing duplicates'
    SV_indicants_no_dups_bam = compilation_prefix + '.no_dups.bam'
    rmdup(SV_indicants_bam, SV_indicants_no_dups_bam)

    # Convert SV indicants bam to sam
    SV_indicants_sam = compilation_prefix + '.no_dups.sam'
    make_sam(SV_indicants_no_dups_bam, SV_indicants_sam)

    # Add mate pairs to SV indicants sam
    print 'adding mate pairs'
    SV_indicants_with_pairs_sam = compilation_prefix + '.with_pairs.sam'
    add_paired_mates(
            SV_indicants_sam, alignment_bam, SV_indicants_with_pairs_sam)

    # Make bam of SV indicants w/mate pairs
    SV_indicants_with_pairs_bam = compilation_prefix + '.with_pairs.bam'
    make_bam(SV_indicants_with_pairs_sam, SV_indicants_with_pairs_bam)

    # Filter low quality reads
    print 'filtering out low quality reads'
    SV_indicants_filtered = compilation_prefix + '.with_pairs.filtered.bam'
    filter_low_qual_read_pairs(
            SV_indicants_with_pairs_bam, SV_indicants_filtered)

    # Sort for velvet assembly
    print 'sorting by name'
    sort_bam_by_name(SV_indicants_filtered)

    return SV_indicants_filtered


def assemble_with_velvet(data_dir, velvet_opts, sv_indicants_bam,
        sample_alignment):

    timestamp = str(datetime.datetime.now())
    contig_number_pattern = re.compile('^NODE_(\d+)_')

    reference_genome = sample_alignment.alignment_group.reference_genome

    contig_files = []
    contig_list = []

    # Write sv_indicants filename and velvet options to file
    assembly_metadata_fn = os.path.join(data_dir, 'metadata.txt')
    with open(assembly_metadata_fn, 'w') as fh:
        assembly_metadata = {
            'sv_indicants_bam': sv_indicants_bam,
            'velvet_opts': velvet_opts
        }
        pickle.dump(assembly_metadata, fh)

    cmd = ' '.join([
            VELVETH_BINARY,
            data_dir,
            str(velvet_opts['velveth']['hash_length'])] +
            ['-' + key + ' ' + str(velvet_opts['velveth'][key])
             for key in velvet_opts['velveth']
             if key not in ['hash_length']] +
            ['-bam', '-shortPaired', sv_indicants_bam])
    print 'velveth cmd:', cmd

    velveth_error_output = os.path.join(data_dir, 'velveth_error_log.txt')
    with open(velveth_error_output, 'w') as error_output_fh:
        subprocess.check_call(cmd, shell=True, executable=settings.BASH_PATH,
                stderr=error_output_fh)

    ins_length = velvet_opts['velvetg'].get('ins_length', None)
    exp_cov = velvet_opts['velvetg'].get('exp_cov', None)
    ins_length_sd = velvet_opts['velvetg'].get('ins_length_sd', None)
    cov_cutoff = velvet_opts['velvetg'].get('cov_cutoff', None)
    min_contig_lgth = velvet_opts['velvetg'].get('min_contig_lgth', None)

    arg_list = [
            VELVETG_BINARY,
            data_dir,
            '-ins_length', str(ins_length),
            '-exp_cov', str(exp_cov),
            '-ins_length_sd', str(ins_length_sd),
            '-cov_cutoff', str(cov_cutoff),
            # '-max_coverage', str(3 * exp_cov),
            '-min_contig_lgth', str(min_contig_lgth),
            '-read_trkg', 'yes']

    cmd = ' '.join(arg_list)
    print 'velvetg cmd:', cmd

    velvetg_error_output = os.path.join(data_dir, 'velveth_error_log.txt')
    with open(velvetg_error_output, 'w') as error_output_fh:
        subprocess.check_call(cmd, shell=True, executable=settings.BASH_PATH,
                stderr=error_output_fh)

    # Collect resulting contigs fasta
    contigs_fasta = os.path.join(data_dir, 'contigs.fa')
    contig_files.append(contigs_fasta)

    records = list(SeqIO.parse(contigs_fasta, 'fasta'))
    digits = len(str(len(records)))
    for seq_record in records:

        if 'N' in seq_record.seq:
            continue

        contig_node_number = int(
                contig_number_pattern.findall(seq_record.description)[0])
        leading_zeros = digits - len(str(contig_node_number))
        contig_label = '%s_%s' % (
                sample_alignment.experiment_sample.label,
                leading_zeros * '0' + str(contig_node_number))

        # Create an insertion model for the contig
        contig = Contig.objects.create(
                label=contig_label,
                parent_reference_genome=reference_genome,
                experiment_sample_to_alignment=(
                        sample_alignment))
        contig_list.append(contig)

        contig.metadata['coverage'] = float(
                seq_record.description.rsplit('_', 1)[1])
        contig.metadata['timestamp'] = timestamp
        contig.metadata['node_number'] = contig_node_number
        contig.metadata['assembly_dir'] = data_dir

        contig.ensure_model_data_dir_exists()
        dataset_path = os.path.join(contig.get_model_data_dir(),
                'fasta.fa')

        with open(dataset_path, 'w') as fh:
            SeqIO.write([seq_record], fh, 'fasta')

        add_dataset_to_entity(
                contig,
                'contig_fasta',
                Dataset.TYPE.REFERENCE_GENOME_FASTA,
                filesystem_location=dataset_path)

    evaluate_contigs(contig_list)
    return contig_files


def evaluate_contigs(contig_list, skip_extracted_read_alignment=False,
        use_read_alignment=True):

    def _length_weighted_coverage(contig):
        return contig.num_bases * contig.coverage

    # Sort contig_list by highest length weighted coverage
    contig_list.sort(key=_length_weighted_coverage, reverse=True)

    # Get placeable contigs using graph-based placement
    placeable_contigs, var_dict_list = graph_contig_placement(contig_list,
            skip_extracted_read_alignment, use_read_alignment)

    # Mark placeable contigs
    for contig in placeable_contigs:
        contig.metadata['is_placeable'] = True

    # Get path for contigs vcf
    contig = contig_list[0]
    sample_alignment = contig.experiment_sample_to_alignment
    ag = contig.experiment_sample_to_alignment.alignment_group
    vcf_path = os.path.join(
            sample_alignment.get_model_data_dir(),
            'de_novo_assembled_contigs.vcf')
    var_dict_vcf_path = os.path.join(
            sample_alignment.get_model_data_dir(),
            'de_novo_assembly_translocations.vcf')

    if var_dict_list:

        dataset_query = sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK)

        if dataset_query:
            dataset_query.delete()

        # Write variant dicts to vcf
        export_var_dict_list_as_vcf(var_dict_list, var_dict_vcf_path,
                contig.experiment_sample_to_alignment)

        # Make dataset for contigs vcf
        add_dataset_to_entity(
                sample_alignment,
                Dataset.TYPE.VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK,
                Dataset.TYPE.VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK,
                var_dict_vcf_path)

    if not placeable_contigs:
        return

    dataset_query = sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.VCF_DE_NOVO_ASSEMBLED_CONTIGS)

    if dataset_query:
        dataset_query.delete()

    # Write contigs to vcf
    export_contig_list_as_vcf(placeable_contigs, vcf_path)

    # Make dataset for contigs vcf
    add_dataset_to_entity(
            sample_alignment,
            Dataset.TYPE.VCF_DE_NOVO_ASSEMBLED_CONTIGS,
            Dataset.TYPE.VCF_DE_NOVO_ASSEMBLED_CONTIGS,
            vcf_path)


@task
def parse_variants_from_vcf(sample_alignment):

    vcf_datasets_to_parse = [
            Dataset.TYPE.VCF_DE_NOVO_ASSEMBLED_CONTIGS,
            Dataset.TYPE.VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK]

    variant_list = []
    for dataset_type in vcf_datasets_to_parse:
        dataset_query = sample_alignment.dataset_set.filter(
            type=Dataset.TYPE.VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK)
        if dataset_query:
            assert dataset_query.count() == 1
            parsed_variants = parse_vcf(
                    dataset_query[0],
                    sample_alignment.alignment_group)

            variant_list.extend(parsed_variants)

    # Add contig origin data to vccds of created variants
    for variant in variant_list:
        vccd_list = variant.variantcallercommondata_set.all()

        for vccd in vccd_list:
            # See if the variant is associated with a contig
            contig_uid_list = vccd.data.get('INFO_contig_uid', False)
            if contig_uid_list:
                contig = Contig.objects.get(uid=contig_uid_list[0])
                contig.variant_caller_common_data = vccd
                contig.save()

        variant.save()
