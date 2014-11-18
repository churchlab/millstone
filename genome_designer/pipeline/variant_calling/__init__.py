"""
Functions for calling Variants.
"""

import os
import shutil
import subprocess
from uuid import uuid4

from celery import task
from django.conf import settings
import vcf

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.s3 import project_files_needed
from pipeline.read_alignment import get_insert_size_mean_and_stdev
from pipeline.variant_calling.lumpy import run_lumpy
from utils.jbrowse_util import add_vcf_track
from utils import uppercase_underscore
from variants.vcf_parser import parse_alignment_group_vcf
from variants.variant_sets import add_variants_to_set_from_bed
from pipeline.variant_effects import run_snpeff

# TODO: These VCF types should be set somewhere else. snpeff_util and
# vcf_parser also use them, but where should they go? settings.py seems
# logical, but it cannot import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF

# Dataset type for results of finding SVs.
VCF_PINDEL_TYPE = Dataset.TYPE.VCF_PINDEL
VCF_DELLY_TYPE = Dataset.TYPE.VCF_DELLY


# Returns a dictionary of common parameters required for all variant callers
# (freebayes, pindel, delly, lumpy).
def get_common_tool_params(alignment_group):
    alignment_type = Dataset.TYPE.BWA_ALIGN
    return {
        'alignment_group': alignment_group,
        'alignment_type': alignment_type,
        'fasta_ref': _get_fasta_ref(alignment_group),
        'output_dir': _create_output_dir(alignment_group),
        'sample_alignments': _find_valid_sample_alignments(
                alignment_group, alignment_type),
    }


def _get_fasta_ref(alignment_group):
    # Grab the reference genome fasta for the alignment.
    return get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()


def _create_output_dir(alignment_group):
    # Prepare a directory to put the output files.
    # We'll put them in
    #     /projects/<project_uid>/alignment_groups/vcf/<variant tool>/
    #     <alignment_type>.vcf
    # We'll save these for now, maybe it's not necessary later.
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    return vcf_dir

def _find_valid_sample_alignments(alignment_group, alignment_type):
    """ Returns a list sample alignment objects for an alignment,
        skipping those that failed. """
    sample_alignment_list = (
            alignment_group.experimentsampletoalignment_set.all())

    # Filter out mis-aligned files.
    # TODO: Should we show in the UI that some alignments failed and are
    # being skipped?
    def _is_successful_alignment(sample_alignment):
        bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)
        return bam_dataset.status == Dataset.STATUS.READY
    sample_alignment_list = [sample_alignment for sample_alignment in
            sample_alignment_list if _is_successful_alignment(sample_alignment)]

    if len(sample_alignment_list) == 0:
        raise Exception('No successful alignments, Freebayes cannot proceed.')

    bam_files = _get_dataset_paths(sample_alignment_list, alignment_type)

    # Keep only valid bam_files
    valid_bam_files = []
    for bam_file in bam_files:
        if bam_file is None:
            continue
        if not os.stat(bam_file).st_size > 0:
            continue
        valid_bam_files.append(bam_file)
    assert len(valid_bam_files) == len(sample_alignment_list), (
            "Expected %d bam files, but found %d" % (
                    len(sample_alignment_list), len(bam_files)))
    return sample_alignment_list


@task
@project_files_needed
def find_variants_with_tool(alignment_group, variant_params):
    """Applies a variant caller to the alignment data contained within
    alignment_group.

    Args:
        alignment_group: AlignmentGroup with all alignments complete.
        variant_params: Triple (tool_name, vcf_dataset_type, tool_function).

    Returns:
        Boolean indicating whether we made it through this entire function.
    """
    # TODO: More informative failure information.
    try:
        common_params = get_common_tool_params(alignment_group)
    except:
        alignment_group.status = AlignmentGroup.STATUS.FAILED
        alignment_group.save(update_fields=['status'])
        return False

    tool_name, vcf_dataset_type, tool_function = variant_params

    # Finding variants means that all the aligning is complete, so now we
    # are VARIANT_CALLING.
    alignment_group.status = AlignmentGroup.STATUS.VARIANT_CALLING
    alignment_group.save()

    # Create subdirectory for this tool
    tool_dir = os.path.join(common_params['output_dir'], tool_name)
    ensure_exists_0775_dir(tool_dir)
    vcf_output_filename = os.path.join(tool_dir,
            uppercase_underscore(common_params['alignment_type']) + '.vcf')

    # Run the tool
    tool_succeeded = tool_function(
            vcf_output_dir= tool_dir,
            vcf_output_filename= vcf_output_filename,
            **common_params)
    if not tool_succeeded:
        return False

    # Add dataset
    # If a Dataset already exists, delete it, might have been a bad run.
    existing_set = Dataset.objects.filter(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    if len(existing_set) > 0:
        existing_set[0].delete()

    vcf_dataset = Dataset.objects.create(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    alignment_group.dataset_set.add(vcf_dataset)

    # Do the following only for freebayes; right now just special if condition
    if tool_name == 'freebayes':
        # For now, automatically run snpeff if a genbank annotation is
        # available.
        # If no annotation, then skip it, and pass the unannotated vcf type.
        if alignment_group.reference_genome.is_annotated():
            run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
            vcf_dataset_type = VCF_ANNOTATED_DATASET_TYPE
        else:
            vcf_dataset_type = VCF_DATASET_TYPE

    sort_vcf(vcf_dataset.get_absolute_location())

    # Tabix index and add the VCF track to Jbrowse
    add_vcf_track(alignment_group.reference_genome, alignment_group,
        vcf_dataset_type)

    # Parse the resulting vcf, grab variant objects
    parse_alignment_group_vcf(alignment_group, vcf_dataset_type)

    flag_variants_from_bed(alignment_group, Dataset.TYPE.BED_CALLABLE_LOCI)

    return True


def flag_variants_from_bed(alignment_group, bed_dataset_type):

    sample_alignments = alignment_group.experimentsampletoalignment_set.all()
    for sample_alignment in sample_alignments:

        # If there is no callable_loci bed, skip the sample alignment.
        # TODO: Make this extensible to other BED files we might have
        callable_loci_bed = get_dataset_with_type(
                entity=sample_alignment,
                type=Dataset.TYPE.BED_CALLABLE_LOCI)

        if not callable_loci_bed: continue

        # need to add sample_alignment and bed_dataset here.
        add_variants_to_set_from_bed(
                sample_alignment= sample_alignment,
                bed_dataset= callable_loci_bed)


def run_freebayes(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, **kwargs):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.

    Returns:
        Boolean, True if successfully made it to the end, else False.
    """
    vcf_dataset_type = VCF_DATASET_TYPE

    bam_files = _get_dataset_paths(sample_alignments, alignment_type)

    # Build up the bam part of the freebayes binary call.
    bam_part = []
    for bam_file in bam_files:
        bam_part.append('--bam')
        bam_part.append(bam_file)

    alignment_group = sample_alignments[0].alignment_group

    if alignment_group.alignment_options['call_as_haploid']:
        alignment_ploidy = 1
    else:
        alignment_ploidy = 2

    # Build the full command and execute it for all bam files at once.
    full_command = (['%s/freebayes/freebayes' %  settings.TOOLS_DIR] + bam_part + [
        '--fasta-reference', fasta_ref,
        '--pvar', '0.001',
        '--ploidy', str(alignment_ploidy),
        '--min-alternate-fraction', '.3',
        '--hwe-priors-off',
        #'--binomial-obs-priors-off',
        '--use-mapping-quality',
        '--min-base-quality', '25',
        '--min-mapping-quality', '30'
    ])

    with open(vcf_output_filename, 'w') as fh:
        subprocess.check_call(full_command, stdout=fh)

    return True # success


def run_pindel(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, **kwargs):
    """Run pindel to find SVs."""
    vcf_dataset_type = VCF_PINDEL_TYPE

    if not os.path.isdir('%s/pindel' % settings.TOOLS_DIR):
        raise Exception('Pindel is not installed. Aborting.')

    bam_files = _get_dataset_paths(sample_alignments, alignment_type)
    samples = [sa.experiment_sample for sa in sample_alignments]
    insert_sizes = [get_insert_size_mean_and_stdev(sa) for sa in
            sample_alignments]

    assert len(bam_files) == len(insert_sizes)

    # Create pindel config file
    pindel_config = os.path.join(vcf_output_dir, 'pindel_config.txt')
    at_least_one_config_line_written = False
    with open(pindel_config, 'w') as fh:
        for bam_file, sample, insert_size in zip(
                bam_files, samples, insert_sizes):

            # Skip bad alignments.
            mean, stdev = insert_size
            if mean == -1:
                continue
            fh.write('%s %s %s\n' % (bam_file, mean, sample.uid))
            at_least_one_config_line_written = True

    if not at_least_one_config_line_written:
        raise Exception
        return False # failure

    # Build the full pindel command.
    pindel_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    subprocess.check_call(['%s/pindel/pindel' % settings.TOOLS_DIR,
        '-f', fasta_ref,
        '-i', pindel_config,
        '-c', 'ALL',
        '-o', pindel_root
    ])

    # convert all different structural variant types to vcf
    subprocess.check_call(['%s/pindel/pindel2vcf' % settings.TOOLS_DIR,
        '-P', pindel_root,
        '-r', fasta_ref,
        '-R', 'name',
        '-d', 'date',
        '-mc', '1',  # just need one read to show 1/1 in vcf
    ])

    postprocess_pindel_vcf(vcf_output_filename)

    return True # success


def run_delly(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, **kwargs):
    """Run delly to find SVs."""
    vcf_dataset_type = VCF_DELLY_TYPE

    if not os.path.isdir('%s/delly' % settings.TOOLS_DIR):
        raise Exception('Delly is not installed. Aborting.')

    delly_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    transformations = ['DEL', 'DUP', 'INV']
    vcf_outputs = map(lambda transformation:
            '%s_%s.vcf' % (delly_root, transformation), transformations)

    # Rename bam files, because Delly uses the name of the file as sample uid.
    # Use cp instead of mv, because other sv callers will be reading from the
    #   original bam file.

    bam_files = _get_dataset_paths(sample_alignments, alignment_type)
    samples = [sa.experiment_sample for sa in sample_alignments]

    new_bam_files = []
    for bam_file, sample in zip(bam_files, samples):
        new_bam_file = os.path.join(
                os.path.dirname(bam_file), sample.uid + '.bam')
        subprocess.check_call(['cp', bam_file, new_bam_file])
        subprocess.check_call(['cp', bam_file + '.bai', new_bam_file + '.bai'])
        new_bam_files.append(new_bam_file)

    # run delly for each type of transformation
    for transformation, vcf_output in zip(transformations, vcf_outputs):

        # not checked_call, because delly errors if it doesn't find any SVs
        subprocess.call(['%s/delly/delly' % settings.TOOLS_DIR,
            '-t', transformation,
            '-o', vcf_output,
            '-g', fasta_ref] + new_bam_files)

    # combine the separate vcfs for each transformation
    vcf_outputs = filter(lambda file: os.path.exists(file), vcf_outputs)
    if vcf_outputs:
        temp_vcf = os.path.join(vcf_output_dir, 'temp_vcf')
        with open(temp_vcf, 'w') as fh:
            subprocess.check_call(['vcf-concat'] + vcf_outputs, stdout=fh)
        with open(vcf_output_filename, 'w') as fh:
            subprocess.check_call(['vcf-sort', temp_vcf], stdout=fh)
        subprocess.check_call(['rm', temp_vcf])
    else:
        # hack: create empty vcf
        subprocess.check_call(['touch', delly_root])
        subprocess.check_call(['%s/pindel/pindel2vcf' % settings.TOOLS_DIR,
            '-p', delly_root,  # TODO does this work?
            '-r', fasta_ref,
            '-R', 'name',
            '-d', 'date'
        ])

    # Delete temporary renamed bam files
    for bam_file in new_bam_files:
        subprocess.check_call(['rm', bam_file])
        subprocess.check_call(['rm', bam_file + '.bai'])

    postprocess_delly_vcf(vcf_output_filename)

    return True # success

# Get paths for each of the dataset files.
def _get_dataset_paths(sample_alignment_list, dataset_type):

    dataset_locations = []

    # These sample alignments should have already
    # been validated in _find_valid_sample_alignments...
    for sample_alignment in sample_alignment_list:
        dataset = get_dataset_with_type(sample_alignment, dataset_type)
        dataset_locations.append(dataset.get_absolute_location())

    return dataset_locations


def _common_postprocess_vcf(vcf_reader):
    # Do postprocessing in common to Pindel and Delly VCFs.
    modified_header_lines = []

    # These properties should be part of VA_DATA, although SV tools will
    #   output at most one property in each row and set Number=1 in the
    #   VCF header line. The easiest way to store these properties in
    #   VA_DATA is just to postprocess the VCF here and change these
    #   header lines to all say Number=A.
    va_properties = ['SVTYPE', 'SVLEN']
    def modify_header(header_line):
        if any([prop in header_line for prop in va_properties]):
            header_line = header_line.replace('Number=1', 'Number=A')
            modified_header_lines.append(header_line)
        return header_line
    vcf_reader._header_lines = map(modify_header, vcf_reader._header_lines)

    # Also add a field for METHOD.
    method_header_line = '##INFO=<ID=METHOD,Number=1,Type=String,' + \
        'Description="Type of approach used to detect SV">\n'
    modified_header_lines.append(method_header_line)
    vcf_reader._header_lines.append(method_header_line)

    # Now update the header lines in vcf_reader.infos map as well.
    parser = vcf.parser._vcf_metadata_parser()
    for header_line in modified_header_lines:
        key, val = parser.read_info(header_line)
        vcf_reader.infos[key] = val


# Postprocess vcfs output by Pindel and Delly, so their information is
#   customized to whatever is needed in Millstone, and the format is
#   the same as that of Freebayes.
def postprocess_pindel_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file))

    _common_postprocess_vcf(vcf_reader)

    # Write the modified VCF to a temp file.
    vcf_writer = vcf.Writer(open(vcf_file + '.tmp', 'a'), vcf_reader)
    for record in vcf_reader:
        if 'SVLEN' not in record.__dict__['INFO']:
            continue  # should not happen

        # pindel uses negative SVLEN for deletions; make them positive
        # always have one entry
        svlen = abs(record.__dict__['INFO']['SVLEN'][0])
        record.__dict__['INFO']['SVLEN'] = [svlen]

        if svlen < 10:  # ignore small variants
            continue

        # update METHOD field
        record.__dict__['INFO']['METHOD'] = 'PINDEL'

        vcf_writer.write_record(record)

    # move temporary file back
    subprocess.check_call(['mv', vcf_file + '.tmp', vcf_file])

def postprocess_delly_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file))

    _common_postprocess_vcf(vcf_reader)

    vcf_writer = vcf.Writer(open(vcf_file + '.tmp', 'a'), vcf_reader)
    for record in vcf_reader:
        record.__dict__['INFO']['METHOD'] = 'DELLY'
        vcf_writer.write_record(record)

    subprocess.check_call(['mv', vcf_file + '.tmp', vcf_file])


# Map from variant caller name to params.
VARIANT_TOOL_PARAMS_MAP = {
    'freebayes': {
        'dataset_type': Dataset.TYPE.VCF_FREEBAYES,
        'runner_fn': run_freebayes
    },
    'pindel': {
        'dataset_type': Dataset.TYPE.VCF_PINDEL,
        'runner_fn': run_pindel
    },
    'delly': {
        'dataset_type': Dataset.TYPE.VCF_DELLY,
        'runner_fn': run_delly
    },
    'lumpy': {
        'dataset_type': Dataset.TYPE.VCF_LUMPY,
        'runner_fn': run_lumpy
    }
}


def get_variant_tool_params():
    """DEPRECATED: Use VARIANT_TOOL_PARAMS_MAP directly.

    Returns a tuple of variant tools params to pass into
    find_variants_with_tool().
    """
    tool_params_list = []
    for tool in ('freebayes', 'pindel', 'delly', 'lumpy'):
        tool_params = VARIANT_TOOL_PARAMS_MAP[tool]
        tool_params_list.append(
                adapt_variant_tool_params_to_tuple(tool, tool_params))
    return tool_params_list


def adapt_variant_tool_params_to_tuple(tool_name, params_dict):
    return (tool_name, params_dict['dataset_type'], params_dict['runner_fn'])


def sort_vcf(input_vcf_filepath):
    """Sorts a vcf file by chromosome and position.

    Overwrites the input.
    """
    temp_vcf = os.path.splitext(input_vcf_filepath)[0] + str(uuid4())[:8]
    assert not os.path.exists(temp_vcf)

    sort_cmd = (
            '(grep ^"#" {original_vcf}; grep -v ^"#" {original_vcf} | '
            'sort -k1,1 -k2,2n) > {sorted_vcf}'
    ).format(
            original_vcf=input_vcf_filepath,
            sorted_vcf=temp_vcf
    )
    subprocess.call(sort_cmd, shell=True)

    shutil.move(temp_vcf, input_vcf_filepath)
