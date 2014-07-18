"""Functions for attempting to de novo align unmapped reads.
"""

import os
import subprocess

from django.conf import settings

from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.models import Dataset


def get_unmapped_reads(sample_alignment, force_rerun=False):
    """Return Dataset for unmapped reads from the sample alignment.
    """
    # First, check if completed dataset already exists.
    unmapped_dataset = get_dataset_with_type(
            sample_alignment, Dataset.TYPE.BWA_UNMAPPED)
    if not force_rerun and unmapped_dataset is not None:
        if (unmapped_dataset.status == Dataset.STATUS.READY and
                os.path.exists(unmapped_dataset.get_absolute_location())):
            return unmapped_dataset
    else:
        unmapped_dataset = Dataset.objects.create(
            label=Dataset.TYPE.BWA_UNMAPPED,
            type=Dataset.TYPE.BWA_UNMAPPED,
            status=Dataset.STATUS.NOT_STARTED)
        sample_alignment.dataset_set.add(unmapped_dataset)

    # If here, we are going to run or re-run the Dataset.
    unmapped_dataset.status = Dataset.STATUS.NOT_STARTED
    unmapped_dataset.save(update_fields=['status'])

    # Get the original bam.
    bam_dataset = get_dataset_with_type(sample_alignment,
            Dataset.TYPE.BWA_ALIGN)
    bam_filename = bam_dataset.get_absolute_location()

    # Allocate a filename for the unmapped reads.
    unmapped_reads_bam_file = os.path.splitext(bam_filename)[0] + '.unmapped.bam'

    try:
        # Start computing.
        unmapped_dataset.status = Dataset.STATUS.COMPUTING
        unmapped_dataset.save(update_fields=['status'])

        # Bits 0x4 set for unmapped segments.
        # We use -f flag to keep only those reads.
        cmd = '{samtools} view -f 0x4 {bam_filename}'.format(
                samtools=settings.SAMTOOLS_BINARY,
                bam_filename=bam_filename)
        with open(unmapped_reads_bam_file, 'w') as output_fh:
           subprocess.check_call(cmd, stdout=output_fh, shell=True)

        # Mark as ready.
        unmapped_dataset.filesystem_location = clean_filesystem_location(
                unmapped_reads_bam_file)
        unmapped_dataset.status = Dataset.STATUS.READY

    except subprocess.CalledProcessError:
        unmapped_dataset.filesystem_location = ''
        unmapped_dataset.status = Dataset.STATUS.FAILED

    unmapped_dataset.save()

    return unmapped_dataset
