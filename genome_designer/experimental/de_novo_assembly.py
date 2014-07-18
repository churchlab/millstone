"""Functions for attempting to de novo align unmapped reads.
"""

import os
import subprocess

from django.conf import settings

from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.models import Dataset
from main.models import get_or_create_derived_bam_dataset


def get_unmapped_reads(sample_alignment, force_rerun=False):
    """Returns Dataset for unmapped reads from the sample alignment.

    Computes it if it doesn't exist already.
    """
    def derivation_fn(sample_alignment, unmapped_reads_dataset):
        # Get the original bam file.
        bam_dataset = get_dataset_with_type(sample_alignment,
                Dataset.TYPE.BWA_ALIGN)
        bam_filename = bam_dataset.get_absolute_location()

        # Allocate a filename for the unmapped reads.
        unmapped_reads_bam_file = (os.path.splitext(bam_filename)[0] +
                '.unmapped.bam')
        unmapped_reads_dataset.filesystem_location = clean_filesystem_location(
                unmapped_reads_bam_file)
        unmapped_reads_dataset.save(update_fields=['filesystem_location'])

        cmd = '{samtools} view -f 0x4 {bam_filename}'.format(
                samtools=settings.SAMTOOLS_BINARY,
                bam_filename=bam_filename)
        with open(unmapped_reads_bam_file, 'w') as output_fh:
           subprocess.check_call(cmd, stdout=output_fh, shell=True)

    return get_or_create_derived_bam_dataset(sample_alignment,
            Dataset.TYPE.BWA_UNMAPPED, derivation_fn, force_rerun=force_rerun)
