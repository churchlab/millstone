"""
Utility functions for working with bam files.
"""

import os
import shutil
import subprocess

from django.conf import settings


def filter_bam_file_by_row(input_bam_path, filter_fn, output_bam_path):
    """Filters rows out of a bam file that don't pass a given filter function.

    This function keeps all header lines.

    Args:
        input_bam_path: Absolute path to input bam file.
        filter_fn: Function applied to each row of the input bam and returns a
            Boolean. If True, keeps the row.
        output_bam_path: Absolute path to the output bam file.
    """
    output_root = os.path.splitext(output_bam_path)[0]
    initial_sam_intermediate = output_root + '.sam'
    filtered_sam_intermediate = output_root + '.filtered.sam'
    final_bam = output_root + '.filtered.bam'

    # Convert to SAM (preserve header with -h option).
    with open(initial_sam_intermediate, 'w') as output_fh:
        p_samtools_view = subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-h', input_bam_path],
                stdout=output_fh)

    # Filter.
    with open(filtered_sam_intermediate, 'w') as output_fh:
        with open(initial_sam_intermediate) as input_fh:
            for line in input_fh:
                # Always write header lines.
                if line[0] == '@':
                    output_fh.write(line)
                    continue

                if filter_fn(line):
                    output_fh.write(line)
                    continue

    # Write final bam.
    with open(final_bam, 'w') as fh:
        p_samtools_view = subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-bS',
                        filtered_sam_intermediate],
                stdout=fh)

    # Move temp file to the original file location.
    shutil.move(final_bam, output_bam_path)

    # Delete intermediate files.
    os.remove(initial_sam_intermediate)
    os.remove(filtered_sam_intermediate)
