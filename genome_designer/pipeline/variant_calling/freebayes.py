"""Wrapper for running Freebayes.
"""

import collections
import errno
import fileinput
import glob
import tempfile
import os
import shutil
import subprocess
import vcf

from celery import task
from django.conf import settings

from main.models import Dataset
from main.model_utils import get_dataset_with_type
from pipeline.read_alignment_util import ensure_bwa_index
from pipeline.variant_calling.common import add_vcf_dataset
from pipeline.variant_calling.common import process_vcf_dataset
from pipeline.variant_calling.common import get_common_tool_params
from pipeline.variant_calling.constants import TOOL_FREEBAYES

from pipeline.variant_effects import run_snpeff
from utils import uppercase_underscore

VCF_AF_HEADER = '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternate allele observation frequency, AO/(RO+AO)">'

def freebayes_regions(ref_genome,
        region_size=settings.FREEBAYES_REGION_SIZE):
    """
    Use bamtools (installed as part of freebayes) to intelligently
    generate regions that will be run in freebayes in parallel.

    ref_genome: the reference genome object
    region_size: how many bases each parallelized region 'chunk' will be
    """
    ref_genome_fasta = get_dataset_with_type(ref_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # ensure that the fasta file has an index
    ensure_bwa_index(ref_genome_fasta)

    ref_genome_faidx = ref_genome_fasta + '.fai'

    regions = []

    with open(ref_genome_faidx) as faidx_fh:
        for line in faidx_fh:
            fields = line.strip().split('\t')
            chr_name, chr_len = fields[:2]
            chr_len = int(chr_len)
            end = 0

        while end < chr_len:
            start = end
            end = start + region_size
            if end > chr_len:
                end = chr_len
            regions.append('{chr_name}:{start}-{end}'.format(
                    chr_name=chr_name,
                    start=start,
                    end=end))
            start = end

    return regions


def run_freebayes(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, region=None, **kwargs):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.

    Returns:
        Boolean, True if successfully made it to the end, else False.
    """
    print 'RUNNING FREEBAYES...'

    bam_files = [
            get_dataset_with_type(sa, alignment_type).get_absolute_location()
            for sa in sample_alignments]

    # Build up the bam part of the freebayes binary call.
    bam_part = []
    for bam_file in bam_files:
        bam_part.append('--bam')
        bam_part.append(bam_file)

    # Determine alignment ploidy (haploid or diploid).
    alignment_group = sample_alignments[0].alignment_group
    if alignment_group.alignment_options['call_as_haploid']:
        alignment_ploidy = 1
    else:
        alignment_ploidy = 2

    other_args_part = [
        '--fasta-reference', fasta_ref,
        '--pvar', '0.001',
        '--ploidy', str(alignment_ploidy),
        '--min-alternate-fraction', '.3',
        '--no-population-priors',
        # '--binomial-obs-priors-off',
        '--use-mapping-quality',
        '--min-base-quality', '25',
        '--min-mapping-quality', '30'
    ]

    if region:
        other_args_part.extend(['--region', region])

    # Build the full command and execute it for all bam files at once.
    full_command = (
            ['%s/freebayes/freebayes' % settings.TOOLS_DIR] +
            bam_part +
            other_args_part)

    print ' '.join(full_command)

    # Run Freebayes.
    with open(vcf_output_filename + '.error', 'w') as error_output_fh:
        with open(vcf_output_filename, 'w') as fh:
            subprocess.check_call(
                    full_command, stdout=fh, stderr=error_output_fh)

    # add the allele frequency FORMAT field to the vcf.
    process_freebayes_region_vcf(vcf_output_filename)

    return True # success

def process_freebayes_region_vcf(vcf_output_filename):
    """
    Processes vcf before region merging.

    IF AO and RO are available for an allele, also add alt allele
    percentages (AF), as percentage of total depth can be a good way to filter
    het/hom calls.
    """

    # store the modified VCF in this temporary file, then move it to overwrite
    # the original file when done adding this field.
    temp_fh = tempfile.NamedTemporaryFile(delete=False)

    with open(vcf_output_filename, 'r') as vcf_input_fh:

        vcf_reader = vcf.Reader(vcf_input_fh)

        # Generate extra header row for AF = AO/(RO+AO).
        vcf_reader._header_lines.append(VCF_AF_HEADER)

        key, val = vcf.parser._vcf_metadata_parser().read_format(VCF_AF_HEADER)
        vcf_reader.formats[key] = val

        # A list of all the FORMAT genotype keys, in order
        format_keys = vcf_reader.formats.keys()

        vcf_writer = vcf.Writer(temp_fh, vcf_reader)

        # Write the old records with the new AF FORMAT field
        for record in vcf_reader:

            # This simply appends ':AF' to the record format field
            record.add_format('AF')

            # check if there are multiple alternate alleles
            multi_alts = len(record.ALT) > 1

            for sample in record.samples:

                # Get alt allele frequencies for each alternate allele.
                try:
                    # TODO: Right now, summing multiple alternate alleles because
                    # we turn arrays into strings in the UI.
                    if multi_alts:
                        total_obs = float(sum(sample['AO']) + sample['RO'])

                        if total_obs > 0:
                            af = sum([float(ao) / total_obs for ao in sample['AO']])

                    # if a single alternate allele:
                    else:
                        total_obs = float(sample['AO'] + sample['RO'])

                        if total_obs > 0:
                            af = float(sample['AO']) / total_obs
                except:
                    af = 0.0

                # new namedtuple with the additional format field
                CallData = collections.namedtuple(
                        'CallData',
                        sample.data._fields+('AF',))

                sample.data = CallData(*sample.data, AF=af)

            vcf_writer.write_record(record)

    # close the writer and move the temp file over the original to replace it
    vcf_writer.close()

    shutil.move(temp_fh.name, vcf_output_filename)
    print 'moved from {} to {}'.format(temp_fh.name, vcf_output_filename)


def merge_freebayes_parallel(alignment_group):
    """
    Merge, sort, and make unique all regional freebayes variant calls after
    parallel execution.

    Returns the Dataset pointing to the merged vcf file. If no freebayes files,
    returns None.
    """
    # First, grab all freebayes parallel vcf files.
    common_params = get_common_tool_params(alignment_group)
    partial_freebayes_vcf_output_dir = os.path.join(
            common_params['output_dir'], 'freebayes')

    # Glob all the parial (region-specific) vcf files.
    # Assert that there is at least one.
    vcf_output_filename_prefix = os.path.join(partial_freebayes_vcf_output_dir,
            uppercase_underscore(common_params['alignment_type']) +
            '.partial.*.vcf')
    vcf_files = glob.glob(vcf_output_filename_prefix)
    if not len(vcf_files):
        return None

    # Generate output filename.
    vcf_ouput_filename_merged = os.path.join(partial_freebayes_vcf_output_dir,
            uppercase_underscore(common_params['alignment_type']) + '.vcf')
    vcf_ouput_filename_merged_fh = open(vcf_ouput_filename_merged, 'w')

    streamsort_cmd = ' '.join([
            settings.VCFSTREAMSORT_BINARY,
            '-w 1000 | ',
            settings.VCFUNIQ_BINARY])

    # create a pipe to write to that will sort all the sub-vcfs
    stream_merge_proc = subprocess.Popen(streamsort_cmd,
            stdin=subprocess.PIPE,
            stdout=vcf_ouput_filename_merged_fh,
            shell=True)

    # concatenate all the vcf files w/ fileinput and keep all but the first
    # header and write to the stream_merge_proc stdin pipe
    header=True
    for line in fileinput.input(vcf_files):
        #https://gist.github.com/waylan/2353749
        try:
            if line.startswith('##'):
                if header:
                    stream_merge_proc.stdin.write(line)
                continue
            elif line.startswith('#'):
                if header:
                   stream_merge_proc.stdin.write(line)
                   header=False
                continue
            stream_merge_proc.stdin.write(line)
        except IOError as e:
            if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
                # Stop loop on "Invalid pipe" or "Invalid argument".
                # No sense in continuing with broken pipe.
                break
            else:
                # Raise any other error.
                raise

    # close the stdin to the stream merge, wait for it to finish,
    # and close the file.
    stream_merge_proc.stdin.close()
    stream_merge_proc.wait()
    vcf_ouput_filename_merged_fh.close()

    vcf_dataset_type = Dataset.TYPE.VCF_FREEBAYES

    # add unannotated vcf dataset first
    vcf_dataset = add_vcf_dataset(alignment_group, vcf_dataset_type,
            vcf_ouput_filename_merged)

    # If genome is annotated then run snpeff now,
    # then update the vcf_output_filename and vcf_dataset_type.
    if alignment_group.reference_genome.is_annotated():

        vcf_ouput_filename_merged_snpeff = run_snpeff(
                alignment_group, TOOL_FREEBAYES)

        vcf_dataset_type = Dataset.TYPE.VCF_FREEBAYES_SNPEFF

        vcf_dataset = add_vcf_dataset(
                alignment_group, vcf_dataset_type,
                vcf_ouput_filename_merged_snpeff)

    # generate variants, process, etc
    process_vcf_dataset(alignment_group, vcf_dataset_type)

    #remove the partial vcfs
    for filename in vcf_files:
        os.remove(filename)

    return vcf_dataset
