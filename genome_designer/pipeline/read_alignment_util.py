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

from django.conf import settings

SAMTOOLS_BINARY = settings.SAMTOOLS_BINARY
TOOLS_DIR = settings.TOOLS_DIR

def ensure_bwa_index(ref_genome_fasta, error_output=None):
    """Creates the reference genome index required by bwa, if it doesn't exist
    already.

    We rely on the convention that the index file location is the fasta
    location with the extension '.bwt' appended to it.
    """
    if not os.path.exists(ref_genome_fasta + '.bwt'):
        build_bwa_index(ref_genome_fasta, error_output)

    # Also build the fasta index.
    if not os.path.exists(ref_genome_fasta + '.fai'):
        subprocess.check_call([
            SAMTOOLS_BINARY,
            'faidx',
            ref_genome_fasta
        ], stderr=error_output)

    ref_genome_dict_path = os.path.splitext(ref_genome_fasta)[0] + '.dict'
    if not os.path.exists(ref_genome_dict_path):
        subprocess.check_call([
            'java', '-Xmx1024M',
            '-jar', '%s/picard/CreateSequenceDictionary.jar' % TOOLS_DIR,
            'R=', ref_genome_fasta,
            'O=', ref_genome_dict_path,
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
