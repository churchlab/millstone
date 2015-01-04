"""Wrapper for running Freebayes.
"""

import errno
import fileinput
import glob
import os
import subprocess

from django.conf import settings

from main.models import Dataset
from main.model_utils import get_dataset_with_type
from pipeline.read_alignment_util import ensure_bwa_index
from pipeline.variant_calling.common import add_vcf_dataset
from pipeline.variant_calling.common import process_vcf_dataset
from pipeline.variant_calling.common import get_common_tool_params

from pipeline.variant_effects import run_snpeff
from utils import uppercase_underscore


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
        '--hwe-priors-off',
        # '--binomial-obs-priors-off',
        '--use-mapping-quality',
        '--min-base-quality', '25',
        '--min-mapping-quality', '30'
    ]

    if region:
        other_args_part.extend(['--region',region])

    # Build the full command and execute it for all bam files at once.
    full_command = (
            ['%s/freebayes/freebayes' % settings.TOOLS_DIR] +
            bam_part +
            other_args_part)

    with open(vcf_output_filename, 'w') as fh:
        subprocess.check_call(full_command, stdout=fh)

    return True # success

def merge_freebayes_parallel(alignment_group):
    """
    Merge, sort, and make unique all regional freebayes variant calls after
    parallel execution.

    Returns the path to the merged vcf file.
    """

    # First, grab all freebayes parallel vcf files.
    common_params = get_common_tool_params(alignment_group)
    tool_dir = os.path.join(common_params['output_dir'], 'freebayes')

    vcfuniq_path = settings.VCFUNIQ_BINARY
    vcfstreamsort_path = settings.VCFSTREAMSORT_BINARY

    vcf_output_filename_prefix = os.path.join(tool_dir,
            uppercase_underscore(common_params['alignment_type']) +
            '.partial.*.vcf')

    vcf_ouput_filename_merged = os.path.join(tool_dir,
            uppercase_underscore(common_params['alignment_type']) + '.vcf')
    vcf_ouput_filename_merged_fh = open(vcf_ouput_filename_merged, 'w')

    streamsort_cmd = ' '.join([
            vcfstreamsort_path, '-w 1000 | ', vcfuniq_path])

    # create a pipe to write to that will sort all the sub-vcfs
    stream_merge_proc = subprocess.Popen(streamsort_cmd,
            stdin=subprocess.PIPE,
            stdout=vcf_ouput_filename_merged_fh,
            shell=True)

    # glob all the vcfs
    vcf_files = glob.glob(vcf_output_filename_prefix)

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
    vcf_annotated_dataset_type = Dataset.TYPE.VCF_FREEBAYES_SNPEFF

    # add unannotated vcf dataset first
    vcf_dataset = add_vcf_dataset(alignment_group, vcf_dataset_type,
            vcf_ouput_filename_merged)

    # If genome is annotated then run snpeff now,
    # then update the vcf_output_filename and vcf_dataset_type.
    if alignment_group.reference_genome.is_annotated():

        vcf_ouput_filename_merged_snpeff = run_snpeff(
                alignment_group, Dataset.TYPE.BWA_ALIGN)
        
        vcf_dataset_type = vcf_annotated_dataset_type

        vcf_dataset = add_vcf_dataset(alignment_group, 
                vcf_dataset_type,
                vcf_ouput_filename_merged_snpeff)

    # generate variants, process, etc
    process_vcf_dataset(alignment_group, vcf_dataset_type)

    #remove the partial vcfs
    for filename in vcf_files:
        os.remove(filename)

    return vcf_dataset


