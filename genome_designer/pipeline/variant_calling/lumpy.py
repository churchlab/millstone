"""Functions for running lumpy and processing output.
"""

import subprocess

from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_split_reads


def run_lumpy(
        fasta_ref, sample_alignments, vcf_output_dir, vcf_output_filename,
        alignment_type, **kwargs):
    """Runs lumpy.
    """
    # NOTE: Only supporting single sample alignment for now. Previously we
    # tried to use lumpy for multiple sample alignments but the machine would
    # run out of memory so we are going to limit functionality to single
    # alignment only for now.
    assert len(sample_alignments) == 1

    # Get the .bam file for the alignment.
    sa = sample_alignments[0]
    bam_dataset = get_dataset_with_type(sa, Dataset.TYPE.BWA_ALIGN)
    bam_filename = bam_dataset.get_absolute_location()

    # Get or create discordant reads.
    bam_disc_dataset = get_discordant_read_pairs(sa)
    bam_disc_file = bam_disc_dataset.get_absolute_location()

    # Get or create split reads.
    bam_sr_dataset = get_split_reads(sa)
    bam_sr_file = bam_sr_dataset.get_absolute_location()

    # Run Lumpy Express.
    subprocess.check_call([
        settings.LUMPY_EXPRESS_BINARY,
        '-B', bam_filename,
        '-S', bam_sr_file,
        '-D', bam_disc_file,
        '-o', vcf_output_filename
    ])

    return True  # success


def filter_lumpy_vcf(original_vcf_path, new_vcf_path):
    """Filters lumpy vcf to get rid of noisy values.

    Args:
        original_vcf_path: Full path to starting vcf.
        new_vcf_path: Path where new vcf will be written.

    NOTE: Now that we are using lumpyexpress, we might not need this anymore.
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
