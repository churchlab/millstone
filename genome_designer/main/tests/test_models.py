"""
Tests for models.py.
"""

import os

from django.test import TestCase

from main.models import auto_generate_short_name
from main.models import clean_filesystem_location
from main.models import Project
from main.models import User
from scripts.import_util import import_reference_genome_from_local_file

import settings

TEST_PROJECT_NAME = 'testModels_project'
TEST_REF_GENOME_NAME = 'mg1655_partial'
TEST_REF_GENOME_PATH =  os.path.join(settings.PWD,
    'test_data/full_vcf_test_set/mg1655_tolC_through_zupT.gb')


class TestModels(TestCase):

    def setUp(self):
        """Override.
        """

        TEST_USERNAME = 'testuser'
        TEST_PASSWORD = 'password'
        TEST_EMAIL = 'test@example.com'

        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)

        self.test_project = Project.objects.create(
            title=TEST_PROJECT_NAME,
            owner=user.get_profile())

        self.test_ref_genome = import_reference_genome_from_local_file(
            self.test_project,
            TEST_REF_GENOME_NAME,
            TEST_REF_GENOME_PATH,
            'genbank')

    def test_clean_filesystem_location(self):
        FAKE_ABS_ROOT = '/root/of/all/evil'
        EXPECTED_CLEAN_URL = 'projects/blah'
        dirty_full_url = os.path.join(FAKE_ABS_ROOT, settings.MEDIA_ROOT,
                EXPECTED_CLEAN_URL)
        self.assertEqual(EXPECTED_CLEAN_URL,
                clean_filesystem_location(dirty_full_url))

    def test_auto_generate_short_name(self):
        LONG_NAME = 'E Coli K12 Substrain MG1655'
        EXPECTED_SHORT_NAME = 'e_coli_k12_s'
        self.assertEqual(EXPECTED_SHORT_NAME,
                auto_generate_short_name(LONG_NAME))

    def testSnpeffOnCreateRefGenome(self):
        """Ensure that Snpeff database is created successfully when creating
           a new reference genome object.
        """

        # check that the genbank file was symlinked
        assert os.path.exists(os.path.join(
                self.test_ref_genome.get_snpeff_directory_path(),
                'genes.gb'))

        # check that the db was made
        assert os.path.exists(os.path.join(
                self.test_ref_genome.get_snpeff_directory_path(),
                'snpEffectPredictor.bin'))
