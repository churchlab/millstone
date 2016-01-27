"""
Tests for genome finishing features
"""
import os
import tempfile

from django.conf import settings
from django.test import TestCase

from pipeline.callable_loci import get_callable_loci


GF_TEST_DIR = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test')

class TestCallableLoci(TestCase):

    def test_basic(self):

        test_bam = os.path.join(GF_TEST_DIR,
                'small_mg1655_data/1kb_ins_del_1000/bwa_align.sorted.withmd.bam')

        tdir = tempfile.mkdtemp(
            prefix='filetest_',
        )
        bed_output_path = os.path.join(tdir, 'callable_loci.bed')

        get_callable_loci(test_bam, bed_output_path)
