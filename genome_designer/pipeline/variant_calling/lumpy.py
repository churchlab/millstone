"""Functions for running lumpy and processing output.
"""

import os
import subprocess

from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_insert_size_mean_and_stdev
from pipeline.read_alignment import get_read_length
from pipeline.read_alignment import get_split_reads


# 1. chromosome 1
# 2. interval 1 start
# 3. interval 1 end
# 4. chromosome 2
# 5. interval 2 start
# 6. interval 2 end
# 7. id
# 8. evidence set score
# 9. strand 1
# 10. strand 2
# 11. type
# 12. id of samples containing evidence for this breakpoint
# 13. strand configurations observed in the evidence set
# 14. point within the two breakpoint with the maximum probability
# 15. segment of each breakpoint that contains 95% of the probability
LUMPY_FIELD_NAMES = [
    'chr_1',
    'ivl_1_start',
    'ivl_1_end',
    'chr_2',
    'ivl_2_start',
    'ivl_2_end',
    'id',
    'evidence_score',
    'strand_1',
    'strand_2',
    'svtype',
    'sample_ids',
    'strand_configs',
    'breakpoint_max',
    'breakpoint_95_reg'
]

LUMPY_VCF_HEADER = "\n".join([
    '##fileformat=VCFv4.0',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=END_CHR,Number=-1,Type=String,Description="End chromosome of the variant described in this record">',
    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
    '##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=SVTYPE,Number=-1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=STRAND_1,Number=1,Type=String,Description="Strand Orientation of SV Start">',
    '##INFO=<ID=STRAND_2,Number=1,Type=String,Description="Strand Orientation of SV End">',
    '##INFO=<ID=METHOD,Number=1,Type=String,Description="SV Caller used to predict">',
    '##INFO=<ID=DP,Number=1,Type=String,Description="combined depth across samples">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=INS,Description="Insertion of novel sequence">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AO,Number=1,Type=Integer,Description="Alternate Allele Observations">'
])


LUMPY_VCF_COL_HEADER_COMMON_FIELDS = [
    '#CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
    'FILTER',
    'INFO',
    'FORMAT'
]

VCF_FORMAT_STR = '\t'.join(['{chr_1}', '{ivl_1}',
        '.', 'N', '<{svtype}>', '{evidence_score}',
        '.', '{info_string}', 'GT:AO', '{genotypes}'])

INFO_FORMAT_STR = (
    'IMPRECISE;'
    'SVTYPE={svtype};'
    'END={ivl_2};'
    'END_CHR={chr_2};'
    'STRAND_1={strand_1};'
    'STRAND_2={strand_2};'
    'SVLEN={svlen};'
    'METHOD=LUMPY;'
    'DP={reads}')



def run_lumpy(fasta_ref, sample_alignments, vcf_output_dir,
              vcf_output_filename, alignment_type, **kwargs):
    """Runs lumpy on the alignments.

    See: https://github.com/arq5x/lumpy-sv
    """
    # TODO, look for these three in kwargs before setting?
    global_lumpy_options = ['-mw', 1, '-tt', 0.0]

    # Split reads should be given a higher weight than discordant pairs.
    shared_sr_options = {
        'back_distance': 20,
        'weight': 2,
        'min_mapping_threshold': 20
    }

    shared_pe_options = {
        'back_distance': 20,
        'weight': 1,
        'min_mapping_threshold': 20,
        'discordant_z': 4
    }

    pe_sample_options = []
    sr_sample_options = []

    # lumpy uses integer sample IDs, so this is a lookup table.
    sample_id_dict = {}
    sample_uid_order = []

    # build up the option strings for each sample
    for i, sa in enumerate(sample_alignments):

        sample_uid = sa.experiment_sample.uid
        sample_uid_order.append(sample_uid)
        sample_id_dict[i] = sample_uid

        bam_pe_dataset = get_discordant_read_pairs(sa)
        bam_pe_file = bam_pe_dataset.get_absolute_location()

        bam_sr_dataset = get_split_reads(sa)
        bam_sr_file = bam_sr_dataset.get_absolute_location()

        ins_size, ins_stdev = get_insert_size_mean_and_stdev(sa)
        histo_dataset = get_dataset_with_type(sa,
                Dataset.TYPE.LUMPY_INSERT_METRICS_HISTOGRAM)
        assert histo_dataset, "Histogram could not be computed."
        histo_file = histo_dataset.get_absolute_location()

        sa_bam_file = get_dataset_with_type(sa,
                Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        read_length = get_read_length(sa_bam_file)

        # Build up paired end part of command.
        if bam_pe_dataset.status != Dataset.STATUS.FAILED and bam_pe_file:
            assert os.path.isfile(bam_pe_file), (
                '{file} is not empty but is not a file!'.format(
                    file=bam_sr_file))
            pe_sample_str = ','.join([
                'bam_file:' + bam_pe_file,
                'histo_file:' + histo_file,
                'mean:' + str(ins_size),
                'stdev:' + str(ins_stdev),
                'read_length:' + str(read_length),
                'min_non_overlap:' + str(read_length),
                'id:' + str(i),
            ] + ['%s:%d' % (k, v) for k, v in shared_pe_options.items()])
            pe_sample_options.extend(['-pe', pe_sample_str])
        else:
            pe_sample_str = ''

        # Build split-read part of data.
        if bam_sr_dataset.status != Dataset.STATUS.FAILED and bam_sr_file:
            assert os.path.isfile(bam_sr_file), (
                '{file} is not empty but is not a file!'.format(
                    file=bam_sr_file))

            sr_sample_str = ','.join([
                'bam_file:' + bam_sr_file,
                'id:' + str(i),
            ] + ['%s:%d' % (k, v) for k, v in shared_sr_options.items()])
        else:
            sr_sample_str = ''

        if bam_pe_dataset.status != Dataset.STATUS.FAILED and bam_pe_file:
            pe_sample_options.extend(['-pe', pe_sample_str])
        if bam_sr_dataset.status != Dataset.STATUS.FAILED and bam_sr_file:
            sr_sample_options.extend(['-sr', sr_sample_str])

    # combine the options and call lumpy
    combined_options = [str(o) for o in (global_lumpy_options +
                                         pe_sample_options + sr_sample_options)]

    lumpy_output = (os.path.splitext(vcf_output_filename)[0] + '_' +
                    sample_uid + '_' + '.txt')

    lumpy_cmd_list = ['%s/lumpy/lumpy' % settings.TOOLS_DIR] + combined_options

    try:
        with open(lumpy_output, 'w') as fh:
            subprocess.check_call(lumpy_cmd_list, stdout=fh)
    except subprocess.CalledProcessError, e:
        print 'Lumpy failed. Command: {command}'.format(
                command=' '.join(lumpy_cmd_list))
        raise e

    return convert_lumpy_output_to_vcf(sample_id_dict, sample_uid_order,
            lumpy_output, vcf_output_filename)


def convert_lumpy_output_to_vcf(sample_id_dict, sample_uid_order, lumpy_output,
        vcf_output_filename):
    lumpy_vcf_col_header = "\t".join(
            LUMPY_VCF_COL_HEADER_COMMON_FIELDS + sample_uid_order)

    # Convert the whole Lumpy output to vcf. We'll filter below.
    unfiltered_vcf = (os.path.splitext(vcf_output_filename)[0] +
            '.unfiltered.vcf')

    with open(lumpy_output, 'r') as lumpy_in:
        with open(unfiltered_vcf, 'w') as lumpy_vcf_out:
            lumpy_vcf_out.write(LUMPY_VCF_HEADER + '\n')
            lumpy_vcf_out.write(lumpy_vcf_col_header + '\n')
            for line in lumpy_in:
                fields = dict(zip(LUMPY_FIELD_NAMES, line.split()))
                lumpy_vcf_out.write(_output_lumpy_line_to_vcf(fields,
                        sample_id_dict, sample_uid_order) + '\n')

    # Filter the results. The reason we first copy the entire output above and
    # then filter is so that we can eventually reuse the vcf filtering method at
    # the expense of a bit of extra writing time. This can be optimized later if
    # it becomes a bottleneck.
    # NOTE: This is the first implementation of vcf filtering an at writing only
    # happens for this lumpy tool. Eventually we want to use with other tools.
    filter_lumpy_vcf(unfiltered_vcf, vcf_output_filename)

    return True # success


def _output_lumpy_line_to_vcf(fields, sample_id_dict, sample_uid_order):
    bp_start, bp_end = [i.split(':') for i in
            fields['breakpoint_max'][4:].split(';')]

    _, fields['ivl_1'] = bp_start
    _, fields['ivl_2'] = bp_end

    # SVLEN is only calculable (easily) if it is an intrachromosal event
    if fields['chr_1'] == fields['chr_2']:
        fields['svlen'] = int(fields['ivl_2']) - int(fields['ivl_1'])

        # negative number if it's a deletion
        if 'DEL' in fields['svtype']:
            fields['svlen'] = -1 * fields['svlen']
    else:
        fields['svlen'] == 0

    # looks like TYPE:DELETION
    # first 3 letters, should be INV, DUP, DEL, INS, etc
    fields['svtype'] = fields['svtype'][5:8]

    # Split up the sample ids
    gt_dict = dict([(sample_uid, '.')
            for sample_uid in sample_id_dict.values()])
    fields['sample_ids'] = fields['sample_ids'][4:].split(';')

    # initialize total read count for this SV
    fields['reads'] = 0

    # make the genotype field strings per sample
    for id_field in fields['sample_ids']:
        sample_id, reads = id_field.split(',')
        fields['reads'] += int(reads)
        # convert lumpy sample ID to our UID
        sample_uid = sample_id_dict[int(sample_id)]
        gt_dict[sample_uid] = ':'.join(['1/1', reads])

    fields['genotypes'] = '\t'.join(
        [gt_dict[uid] for uid in sample_uid_order])

    fields['info_string'] = INFO_FORMAT_STR.format(
        **fields)

    return VCF_FORMAT_STR.format(**fields)


def filter_lumpy_vcf(original_vcf_path, new_vcf_path):
    """Filters lumpy vcf to get rid of noisy values.

    Args:
        original_vcf_path: Full path to starting vcf.
        new_vcf_path: Path where new vcf will be written.
    """
    with open(original_vcf_path) as orig_vcf_fh:
        with open(new_vcf_path, 'w') as new_vcf_fh:
            vcf_reader = vcf.Reader(orig_vcf_fh)
            vcf_writer = vcf.Writer(new_vcf_fh, vcf_reader)
            for record in vcf_reader:
                # If record fails any filter, continue to next record without
                # writing.
                if int(record.INFO['DP']) < 10:
                    continue
                vcf_writer.write_record(record)
