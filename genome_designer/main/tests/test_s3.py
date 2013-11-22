"""
Tests for s3.py.
"""

from django.test import TestCase
from django.conf import settings
from django.test.utils import override_settings
import logging
import tempfile
import glob
import os
from main.s3 import *

class TestS3(TestCase):
    def setUp(self):
        if settings.TEST_S3 and settings.S3_DRY_RUN:
            logging.warning("Overriding settings.S3_DRY_RUN to perform S3 tests.")

    def test_calc_file_md5(self):
        t = tempfile.NamedTemporaryFile()
        t.write("A" * 1024 * 1024 * 20)  # 20MBytes

        md5 = calc_file_md5(t.name)
        self.assertEqual(md5, "2047e9836570e1cd4eaf82b8c67fa5a1")

    @override_settings(S3_DRY_RUN=False)
    def test_s3_put_get_string(self):
        if settings.TEST_S3:
            test_string = "blah"
            key = "__tests__/s3_put_string.test"
            s3_put_string(key, test_string)
            self.assertEqual(test_string, s3_get_string(key))

    @override_settings(S3_DRY_RUN=False)
    def test_s3_put_directory(self):
        if settings.TEST_S3:
            s3_test_directory = "__tests__/test_data/fake_genome_and_reads"
            fake_genome_and_reads_dir = os.path.join(settings.PWD,
                'test_data/fake_genome_and_reads')
            self.assertTrue(os.path.exists(fake_genome_and_reads_dir))
            s3_put_directory(s3_test_directory, fake_genome_and_reads_dir)

            files = glob.glob(fake_genome_and_reads_dir + "/**/*")
            for f in files:
                key = aws_bucket.get_key(s3_test_directory + "/" + os.path.relpath(f,
                    fake_genome_and_reads_dir))
                self.assertIsNotNone(key)
                self.assertEqual(calc_file_md5(f), key.etag.strip("\""))

            for key in aws_bucket.list(s3_test_directory):
                key.delete()

    @override_settings(S3_DRY_RUN=False)
    def test_s3_get_directory(self):
        if settings.TEST_S3:
            s3_test_directory = "__tests__/test_data/fake_genome_and_reads"
            fake_genome_and_reads_dir = os.path.join(settings.PWD,
                'test_data/fake_genome_and_reads')
            self.assertTrue(os.path.exists(fake_genome_and_reads_dir))
            s3_put_directory(s3_test_directory, fake_genome_and_reads_dir)

            out_dir = tempfile.mkdtemp()
            s3_get_directory(s3_test_directory, out_dir)

            files = glob.glob(fake_genome_and_reads_dir + "/**/*")
            for f in files:
                out_file = os.path.join(out_dir, os.path.relpath(f,
                    fake_genome_and_reads_dir))
                self.assertTrue(os.path.isfile(out_file))
                self.assertEqual(calc_file_md5(f), calc_file_md5(out_file))

            for key in aws_bucket.list(s3_test_directory):
                key.delete()

    @override_settings(S3_DRY_RUN=False)
    def tearDown(self):
        if settings.TEST_S3:
            for key in aws_bucket.list("__tests__"):
                key.delete()
