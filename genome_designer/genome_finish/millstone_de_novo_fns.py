import subprocess
import os
import pickle

from Bio import SeqIO
from django.conf import settings
import numpy as np
import pysam

from genome_finish import __path__ as gf_path_list
from genome_finish.insertion_placement_read_trkg import extract_left_and_right_clipped_read_dicts
from main.models import Dataset
from main.models import Variant
from main.models import VariantSet
from settings import SAMTOOLS_BINARY
from settings import BASH_PATH
from utils.bam_utils import clipping_stats
from variants.variant_sets import update_variant_in_set_memberships

GENOME_FINISH_PATH = gf_path_list[0]
VELVETH_BINARY = settings.TOOLS_DIR + '/velvet/velveth'
VELVETG_BINARY = settings.TOOLS_DIR + '/velvet/velvetg'


def get_altalign_reads(input_bam_path, output_bam_path, xs_threshold=None):
    input_af = pysam.AlignmentFile(input_bam_path, 'rb')
    output_af = pysam.AlignmentFile(output_bam_path, 'wb',
                template=input_af)

    for read in input_af:
        if read.has_tag('XS') and read.has_tag('AS'):
            if read.get_tag('AS') <= read.get_tag('XS'):
                output_af.write(read)

    output_af.close()
    input_af.close()


def get_clipped_reads(bam_filename, output_filename, clipping_threshold=None):
    """Creates bam of reads that have more than clipping_threshold bases
    of clipping.  If no clipping_threshold specified, clipping stats
    for the alignment are calculated and the clipping_threshold is set
    to the mean + one stddev of the per read clipping of a sample of
    10000 reads
    """
    if clipping_threshold is None:
        stats = clipping_stats(bam_filename, sample_size=10000)
        clipping_threshold = int(stats['mean'] + stats['std'])

    cmd = ' | '.join([
            # -F 0x100 option filters out secondary alignments
            '{samtools} view -h -F 0x100 {bam_filename}',
            '{extract_clipped_script} -i stdin -t {clipping_threshold}',
            '{samtools} view -Sb -']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    clipping_threshold=clipping_threshold,
                    extract_clipped_script=os.path.join(
                            GENOME_FINISH_PATH,
                            'extractClippedReads.py'))

    with open(output_filename, 'w') as fh:
        subprocess.check_call(cmd, stdout=fh, shell=True, executable=BASH_PATH)

    # sort the split reads, overwrite the old file
    subprocess.check_call(
            [SAMTOOLS_BINARY, 'sort', output_filename] +
            [os.path.splitext(output_filename)[0]])


def get_piled_reads(input_bam_path, output_bam_path,
        clipping_threshold=None):
    """Creates bam of reads that have more than clipping_threshold bases
    of clipping and are stacked 3 standard deviations higher than the average
    pileup of clipped reads.  If no clipping_threshold specified, clipping
    stats for the alignment are calculated and the clipping_threshold is set
    to the mean + one stddev of the per read clipping of a sample of
    10000 reads.
    """
    if clipping_threshold is None:
        stats = clipping_stats(input_bam_path, sample_size=10000)
        clipping_threshold = int(stats['mean'] + stats['std'])

    input_af = pysam.AlignmentFile(input_bam_path, 'rb')
    output_af = pysam.AlignmentFile(output_bam_path, 'wb',
            template=input_af)
    lr_clipped = extract_left_and_right_clipped_read_dicts(
            input_af,
            clipping_threshold=clipping_threshold)
    input_af.close()

    for clipped_dict in [
            lr_clipped['left_clipped'], lr_clipped['right_clipped']]:
        stack_counts = map(len, clipped_dict.values())
        mean_stacking = np.mean(stack_counts)
        std_stacking = np.std(stack_counts)
        stacking_cutoff = mean_stacking + 3 * std_stacking
        for read_list in clipped_dict.values():
            if len(read_list) > stacking_cutoff:
                for read in read_list:
                    output_af.write(read)
    output_af.close()


def get_clipped_reads_smart(input_bam_path, output_bam_path,
        clipping_threshold=8, phred_encoding=None):

    """Gets reads not overlapping their adaptor with a terminal
    segment of clipping with average phred scores above the cutoff
    """

    phred_encoding_to_shift = {
        'Illumina 1.5': 31,
        'Sanger / Illumina 1.9': 0
    }

    CLIPPED_AVG_PHRED_CUTOFF = 20
    if (phred_encoding is not None and
        phred_encoding in phred_encoding_to_shift):

        CLIPPED_AVG_PHRED_CUTOFF += phred_encoding_to_shift[phred_encoding]

    SOFT_CLIP = 4
    HARD_CLIP = 5
    CLIP = [SOFT_CLIP, HARD_CLIP]

    input_af = pysam.AlignmentFile(input_bam_path, 'rb')
    output_af = pysam.AlignmentFile(output_bam_path, 'wb',
                template=input_af)

    for read in input_af:
        # If no cigartuples, i.e. unmapped, continue
        if read.cigartuples is None:
            continue

        if read.is_secondary or read.is_supplementary:
            continue

        # TODO: Account for template length
        # adapter_overlap = max(read.template_length - query_alignment_length, 0)

        # Determine left and right clipped counts
        left_clipping = (read.cigartuples[0][1]
                if read.cigartuples[0][0] in CLIP else 0)
        right_clipping = (read.cigartuples[-1][1]
                if read.cigartuples[-1][0] in CLIP else 0)

        # Write reads to file if clipped bases have average phred score
        # above cutoff
        if left_clipping > clipping_threshold:
            clipped_phred_scores = read.query_qualities[:left_clipping]
            if np.mean(clipped_phred_scores) > CLIPPED_AVG_PHRED_CUTOFF:
                output_af.write(read)
                continue
        if right_clipping > clipping_threshold:
            clipped_phred_scores = read.query_qualities[-right_clipping:]
            if np.mean(clipped_phred_scores) > CLIPPED_AVG_PHRED_CUTOFF:
                output_af.write(read)
                continue

    output_af.close()
    input_af.close()


def get_match_counts(bam_filename):
    cmd = ' | '.join([
            '{samtools} view -h {bam_filename}',
            '{assess_alignment_script} -i stdin']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    assess_alignment_script=os.path.join(
                        GENOME_FINISH_PATH,
                        'assess_alignment.py'))

    output = subprocess.check_output(cmd, shell=True, executable=BASH_PATH)
    return pickle.loads(output)


def get_insertion_location_from_bam(bam_path):
    """Attempts to find the insertion locations of a contig
    using the alignment of that contig to the reference
    """

    BAM_CSOFT_CLIP = 4
    BAM_CHARD_CLIP = 5
    clip_codes = [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]

    samfile = pysam.AlignmentFile(bam_path)
    contig_length = 0

    left_read = None
    right_read = None
    for read in samfile:
        contig_length = max(contig_length, read.query_length)
        if read.cigartuples[0][0] not in clip_codes:
            left_read = read
        elif read.cigartuples[-1][0] not in clip_codes:
            right_read = read

    if left_read is None:
        error_string = 'Expected homology on left end of contig not found'
    elif right_read is None:
        error_string = 'Expected homology on right end of contig not found'
    else:
        left_start_pos = left_read.reference_start
        right_end_pos = right_read.reference_end
        chromosome_seqrecord_id = samfile.getrname(left_read.reference_id)

        if left_start_pos > right_end_pos:
            error_string = ('Contig alignment to reference suggests a ' +
                    'structural variant more complex than a simple insertion')
        elif contig_length < right_end_pos - left_start_pos:
            error_string = ('Contig not long enough to span putative ' +
                    'insertion region, suggesting insertion is longer than ' +
                    'assembled contig')
        else:
            locations_dict = {
                    'chromosome_seqrecord_id': chromosome_seqrecord_id,
                    'left_end': left_start_pos,
                    'right_end': right_end_pos,
            }
            return locations_dict

    return {'error_string': error_string}


def make_sliced_fasta(fasta_path, seqrecord_id, left_index, right_index,
        output_fasta_path=None):
    with open(fasta_path, 'r') as fh:
        seqrecord = SeqIO.parse(fh, 'fasta').next()
        while(seqrecord.id != seqrecord_id):
            seqrecord = SeqIO.parse(fh, 'fasta').next()

    seqrecord = seqrecord[left_index:right_index]

    if output_fasta_path:
        fasta_prefix, fasta_ext = os.path.splitext(output_fasta_path)
        sliced_fasta_path = output_fasta_path
    else:
        fasta_prefix, fasta_ext = os.path.splitext(fasta_path)
        sliced_fasta_path = (fasta_prefix + '.split_' + str(left_index) +
                '_to_' + str(right_index) + fasta_ext)
        if os.path.exists(sliced_fasta_path):
            raise Exception("Sliced fasta already exists with path:",
                    sliced_fasta_path)

    with open(sliced_fasta_path, 'w') as fh:
        SeqIO.write(seqrecord, fh, 'fasta')

    return sliced_fasta_path


def get_unmapped_reads(bam_filename, output_filename, avg_phred_cutoff=None):

    if avg_phred_cutoff is not None:
        intermediate_filename = '_unfiltered'.join(
                os.path.splitext(output_filename))
    else:
        intermediate_filename = output_filename

    cmd = '{samtools} view -h -b -f 0x4 {bam_filename}'.format(
            samtools=SAMTOOLS_BINARY,
            bam_filename=bam_filename)
    with open(intermediate_filename, 'w') as output_fh:
        subprocess.call(
                cmd, stdout=output_fh, shell=True, executable=BASH_PATH)

    if avg_phred_cutoff is not None:
        filter_low_qual_read_pairs(intermediate_filename, output_filename,
                avg_phred_cutoff)


def get_split_reads(bam_filename, output_filename):
    """Isolate split reads from a sample alignment.

    This uses a python script supplied with Lumppy, that is run as a
    separate process.

    NOTE THAT THIS SCRIPT ONLY WORKS WITH BWA MEM.
    """

    # Use lumpy bwa-mem split read script to pull out split reads.
    filter_split_reads = ' | '.join([
            # -F 0x100 option filters out secondary alignments
            '{samtools} view -h -F 0x100 {bam_filename}',
            'python {lumpy_bwa_mem_sr_script} -i stdin',
            '{samtools} view -Sb -']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    lumpy_bwa_mem_sr_script=
                            settings.LUMPY_EXTRACT_SPLIT_READS_BWA_MEM)

    try:
        with open(output_filename, 'w') as fh:
            subprocess.check_call(
                    filter_split_reads, stdout=fh, shell=True,
                    executable=BASH_PATH)

        # sort the split reads, overwrite the old file
        subprocess.check_call(
                [SAMTOOLS_BINARY, 'sort', output_filename,
                 os.path.splitext(output_filename)[0]])

    except subprocess.CalledProcessError:
        raise Exception('Exception caught in split reads generator, ' +
                        'perhaps due to no split reads')


def _parse_sam_line(line):
    parts = line.split()
    return {
        'read_id': parts[0],
        'flags': parts[1]
    }


def add_paired_mates(input_sam_path, source_bam_filename, output_sam_path):

    sam_file = pysam.AlignmentFile(input_sam_path)
    input_qnames_to_read = {}
    for read in sam_file:
        input_qnames_to_read[read.qname] = True
    sam_file.close()

    original_alignmentfile = pysam.AlignmentFile(source_bam_filename, "rb")
    output_alignmentfile = pysam.AlignmentFile(
            output_sam_path, "wh", template=original_alignmentfile)
    for read in original_alignmentfile:
        if input_qnames_to_read.get(read.qname, False):
            if not read.is_secondary and not read.is_supplementary:
                output_alignmentfile.write(read)
    output_alignmentfile.close()
    original_alignmentfile.close()


def filter_out_unpaired_reads(input_bam_path, output_bam_path):

    input_af = pysam.AlignmentFile(input_bam_path, 'rb')

    # Build qname -> flag list dictionary
    read_flags = {}
    for read in input_af:
        if read.qname not in read_flags:
            read_flags[read.qname] = [read.flag]
        else:
            read_flags[read.qname].append(read.flag)

    # Build qname -> is_paired dictionary
    reads_with_pairs = {}
    not_primary_alignment_flag = 256
    supplementary_alignment_flag = 2048
    for qname, flags in read_flags.items():
        primary_count = 0
        for f in flags:
            if (not (f & not_primary_alignment_flag) and
                    not (f & supplementary_alignment_flag)):
                primary_count += 1
        if primary_count == 2:
            reads_with_pairs[qname] = True

    # Write reads in input to output if not in bad_quality_names
    output_af = pysam.AlignmentFile(output_bam_path, "wb",
            template=input_af)
    input_af.reset()
    for read in input_af:
        if read.qname in reads_with_pairs:
            output_af.write(read)
    output_af.close()
    input_af.close()


def filter_low_qual_read_pairs(input_bam_path, output_bam_path,
        avg_phred_cutoff=20):
    """Filters out reads with average phred scores below cutoff
    """

    # Put qnames with average phred scores below the cutoff into dictionary
    bad_quality_qnames = {}
    input_af = pysam.AlignmentFile(input_bam_path, "rb")
    for read in input_af:
        avg_phred = np.mean(read.query_qualities)
        if avg_phred < avg_phred_cutoff:
            bad_quality_qnames[read.qname] = True

    # Write reads in input to output if not in bad_quality_names
    output_af = pysam.AlignmentFile(output_bam_path, "wb",
            template=input_af)
    input_af.reset()
    for read in input_af:
        if not bad_quality_qnames.get(read.qname, False):
            output_af.write(read)
    output_af.close()
    input_af.close()


def run_velvet(reads, output_dir, opt_dict=None):
    default_opts = {
        'velveth': {
            'hash_length': 21
        },
    }

    if not opt_dict:
        opt_dict = default_opts

    if 'hash_length' not in opt_dict['velveth']:
        raise Exception('If passing an option dictionary, an output_dir_name key \
        must exist in the velveth sub-dictionary')

    cmd = ' '.join([
            VELVETH_BINARY,
            output_dir,
            str(opt_dict['velveth']['hash_length'])] +
            ['-' + key + ' ' + str(opt_dict['velveth'][key])
             for key in opt_dict['velveth']
             if key not in ['hash_length']] +
            ['-bam', reads])
    subprocess.check_call(cmd, shell=True, executable=BASH_PATH)

    cmd = ' '.join([
            VELVETG_BINARY,
            output_dir] +
            ['-' + key + ' ' + str(opt_dict['velvetg'][key])
             for key in opt_dict['velvetg']])
    subprocess.check_call(cmd, shell=True, executable=BASH_PATH)


def create_de_novo_variants_set(alignment_group, variant_set_label,
        callers_to_not_include=['COV_DETECT_DELETION']):

    ref_genome = alignment_group.reference_genome

    # Get de novo variants
    de_novo_variants = []
    for variant in Variant.objects.filter(
            reference_genome=ref_genome):
        for vccd in variant.variantcallercommondata_set.all():
            if vccd.data.get('INFO_METHOD', None) == 'DE_NOVO_ASSEMBLY' and (
                    vccd.data.get('INFO_CALLER', None) not in (
                            callers_to_not_include)):
                de_novo_variants.append(variant)
                continue

    variant_set = VariantSet.objects.create(
            reference_genome=ref_genome,
            label='de_novo_variants')


    update_variant_in_set_memberships(
        ref_genome,
        [variant.uid for variant in de_novo_variants],
        'add',
        variant_set.uid)

    return variant_set


def get_coverage_stats(sample_alignment):
    """Returns a dictionary with chromosome seqrecord_ids as keys and
    subdictionaries as values.

    Each subdictionary has three keys: length, mean, and std which hold the
    particular chromosome's length, mean read coverage, and standard
    deviation of read coverage
    """

    maybe_chrom_cov_dict = sample_alignment.data.get('chrom_cov_dict', None)
    if maybe_chrom_cov_dict is not None:
        return maybe_chrom_cov_dict

    bam_path = sample_alignment.dataset_set.get(type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    alignment_af = pysam.AlignmentFile(bam_path)
    chrom_list = alignment_af.references
    chrom_lens = alignment_af.lengths

    c_starts = [0]*len(chrom_list)
    c_ends = chrom_lens

    chrom_cov_lists = []
    for chrom, c_start, c_end in zip(chrom_list, c_starts, c_ends):

        chrom_cov_lists.append([])
        cov_list = chrom_cov_lists[-1]
        for pileup_col in alignment_af.pileup(chrom,
                start=c_start, end=c_end, truncate=True):

            depth = pileup_col.nsegments
            cov_list.append(depth)
    alignment_af.close()

    sub_dict_tup_list = zip(
            chrom_lens,
            map(np.mean, chrom_cov_lists),
            map(np.std, chrom_cov_lists))

    sub_dict_list = map(
            lambda tup: dict(zip(['length', 'mean', 'std'], tup)),
            sub_dict_tup_list)

    chrom_cov_dict = dict(zip(chrom_list, sub_dict_list))

    sample_alignment.data['chrom_cov_dict'] = chrom_cov_dict
    sample_alignment.save()
    return chrom_cov_dict


def get_avg_genome_coverage(sample_alignment):
    """Returns a float which is the average genome coverage, calculated as
    the average length-weighted read coverage over all chromosomes
    """

    coverage_stats = get_coverage_stats(sample_alignment)

    len_weighted_coverage = 0
    total_len = 0

    for sub_dict in coverage_stats.values():
        length = sub_dict['length']
        avg_coverage = sub_dict['mean']

        len_weighted_coverage += length * avg_coverage
        total_len += length

    return float(len_weighted_coverage) / total_len
