"""
Methods related to aligning reads.

This module was created because of a circular import issue with celery.
  ...
  File "/home/glebk/Projects/churchlab/genome-designer-v2/genome_designer/pipeline/read_alignment.py", line 9, in <module>
    from main.models import clean_filesystem_location
  File "/home/glebk/Projects/churchlab/genome-designer-v2/genome_designer/main/__init__.py", line 1, in <module>
    import signals
  File "/home/glebk/Projects/churchlab/genome-designer-v2/genome_designer/main/signals.py", line 20, in <module>
    from pipeline.read_alignment import ensure_bwa_index
ImportError: cannot import name ensure_bwa_index
"""

import os
import subprocess

from utils.bam_utils import filter_bam_file_by_row

from django.conf import settings

SAMTOOLS_BINARY = settings.SAMTOOLS_BINARY
TOOLS_DIR = settings.TOOLS_DIR


def has_bwa_index(ref_genome_fasta):
    return os.path.exists(ref_genome_fasta + '.bwt')


def ensure_bwa_index(ref_genome_fasta, error_output=None):
    """Creates the reference genome index required by bwa, if it doesn't exist
    already.

    We rely on the convention that the index file location is the fasta
    location with the extension '.bwt' appended to it.
    """
    if not has_bwa_index(ref_genome_fasta):
        build_bwa_index(ref_genome_fasta, error_output)

    # Also build the fasta index.
    if not os.path.exists(ref_genome_fasta + '.fai'):
        subprocess.check_call([
            SAMTOOLS_BINARY,
            'faidx',
            ref_genome_fasta
        ], stderr=error_output)


def build_bwa_index(ref_genome_fasta, error_output=None):
    """Calls the command that builds the bwa index required for alignment.

    This creates a file in the same directory as ref_genome_fasta, appending
    the extension '.bwt' to the name of the fasta.
    """
    subprocess.check_call([
        '%s/bwa/bwa' % TOOLS_DIR,
        'index',
        '-a',
        'is',
        ref_genome_fasta
    ], stderr=error_output)


def index_bam_file(bam_file, error_output=None):
    subprocess.check_call([
            SAMTOOLS_BINARY,
            'index',
            bam_file,
            ], stderr=error_output)

def extract_split_reads(bam_filename, bam_split_filename):
    """
    Isolate split reads from a bam file.
    This uses a python script supplied with Lumpy that is run as a
    separate process.

    This is an internal function that works directly with files, and
    is called separately by both SV calling and read ref alignment.

    NOTE THAT THIS SCRIPT ONLY WORKS WITH BWA MEM.
    """
    assert os.path.exists(bam_filename), "BAM file '%s' is missing." % (
            bam_filename)

    # Use lumpy bwa-mem split read script to pull out split reads.
    filter_split_reads = ' | '.join([
            '{samtools} view -h {bam_filename}',
            'python {lumpy_bwa_mem_sr_script} -i stdin',
            '{samtools} view -Sb -']).format(
                    samtools=settings.SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    lumpy_bwa_mem_sr_script=
                            settings.LUMPY_EXTRACT_SPLIT_READS_BWA_MEM)

    with open(bam_split_filename, 'w') as fh:
        subprocess.check_call(filter_split_reads,
                stdout=fh,
                shell=True,
                executable=settings.BASH_PATH)

    # sort the split reads, overwrite the old file
    subprocess.check_call([settings.SAMTOOLS_BINARY, 'sort',
            bam_split_filename,
            os.path.splitext(bam_split_filename)[0]])

    _filter_out_interchromosome_reads(bam_split_filename)

def extract_discordant_read_pairs(bam_filename, bam_discordant_filename):
    """Isolate discordant pairs of reads from a sample alignment.
    """
    # Use bam read alignment flags to pull out discordant pairs only
    filter_discordant = ' | '.join([
            '{samtools} view -u -F 0x0002 {bam_filename} ',
            '{samtools} view -u -F 0x0100 - ',
            '{samtools} view -u -F 0x0004 - ',
            '{samtools} view -u -F 0x0008 - ',
            '{samtools} view -b -F 0x0400 - ']).format(
                    samtools=settings.SAMTOOLS_BINARY,
                    bam_filename=bam_filename)

    with open(bam_discordant_filename, 'w') as fh:
        subprocess.check_call(filter_discordant,
                stdout=fh, shell=True, executable=settings.BASH_PATH)

    # sort the discordant reads, overwrite the old file
    subprocess.check_call([settings.SAMTOOLS_BINARY, 'sort', 
            bam_discordant_filename,
            os.path.splitext(bam_discordant_filename)[0]])

    _filter_out_interchromosome_reads(bam_discordant_filename)


def _filter_out_interchromosome_reads(bam_filename, overwrite_input=True):
    """Filters out read pairs which lie on different chromosomes.

    Args:
        bam_filename: Path to bam file.
        overwrite_input: If True, overwrite the input file.
    """
    def is_rnext_same(line):
        parts = line.split('\t')
        rnext_col = parts[6]
        return rnext_col == '='

    if overwrite_input:
        output_bam_path = bam_filename
    else:
        output_bam_path = os.path.splitext(bam_filename)[0] + '.nointerchrom.bam'

    filter_bam_file_by_row(bam_filename, is_rnext_same, output_bam_path)

