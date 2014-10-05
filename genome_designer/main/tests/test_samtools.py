"""
Tests for samtools.
"""

import os
import subprocess

from django.conf import settings
from django.test import TestCase


alignment_unsorted_sam = os.path.join(settings.PWD,
    'test_data/genome_finish_test/bwa_align.alignment.unsorted.sam')
alignment_unsorted_bam = os.path.join(settings.PWD,
    'test_data/genome_finish_test/bwa_align.alignment.unsorted.bam')


class TestSamtools(TestCase):

    def test_samtools(self):

        cmd = "{samtools} view -b -S {alignment_sam} > {alignment_bam}".format(
                samtools=settings.SAMTOOLS_BINARY,
                alignment_sam=alignment_unsorted_sam,
                alignment_bam=alignment_unsorted_bam)

        # print "test_cmd:", cmd
        subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

        bam_file_empty = os.stat(alignment_unsorted_bam).st_size == 0

        assert not bam_file_empty

        os.remove(alignment_unsorted_bam)
