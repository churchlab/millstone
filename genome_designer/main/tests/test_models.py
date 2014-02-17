"""
Tests for models.py.
"""

import os

from django.conf import settings
from django.test import TestCase

from main.models import auto_generate_short_name
from main.models import clean_filesystem_location
from main.models import AlignmentGroup
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import User
from main.models import Dataset
from scripts.import_util import import_reference_genome_from_local_file
import subprocess


TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'
TEST_PROJECT_NAME = 'testModels_project'
TEST_REF_GENOME_NAME = 'mg1655_partial'
TEST_REF_GENOME_PATH = os.path.join(settings.PWD,
    'test_data/full_vcf_test_set/mg1655_tolC_through_zupT.gb')


class TestModels(TestCase):

    def setUp(self):
        """Override.
        """
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


    def testDataSetCompression(self):
        """Make sure data set compression behaves correctly.
        """
        dataset = Dataset.objects.create(
                label='test_dataset', 
                type=Dataset.TYPE.FASTQ1)

        GZIPPED_FASTQ_FILEPATH = os.path.join(settings.PWD, 'test_data',
                'compressed_fastq', 'sample0.simLibrary.1.fq.gz')

        dataset.filesystem_location = clean_filesystem_location(
                    GZIPPED_FASTQ_FILEPATH)

        assert dataset.is_compressed()

        process = subprocess.Popen(
                ('head '+dataset.wrap_if_compressed()+' | wc -l'), 
                shell=True, executable=settings.BASH_PATH, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

        wc_output, errmsg = process.communicate()
        rc = process.returncode

        assert rc == 0, (
                "Compression process returned non-zero exit status: %s" % (
                        errmsg))

        assert int(wc_output) == 10, (
                "Compression failed: %s" % (errmsg))


class TestAlignmentGroup(TestCase):

    def test_get_samples(self):
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.test_project = Project.objects.create(
                title=TEST_PROJECT_NAME,
                owner=user.get_profile())
        self.test_ref_genome = ReferenceGenome.objects.create(
                project=self.test_project,
                label='blah',
                num_chromosomes=1,
                num_bases=1000)
        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.test_ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

        # Create a bunch of samples and relate them.
        for sample_idx in range(10):
            sample = ExperimentSample.objects.create(
                    uid=str(sample_idx),
                    project=self.test_project,
                    label='some label'
            )
            ExperimentSampleToAlignment.objects.create(
                    alignment_group=alignment_group,
                    experiment_sample=sample)

        # Test the method.
        samples = alignment_group.get_samples()
        sample_uid_set = set([sample.uid for sample in samples])
        self.assertEqual(sample_uid_set,
                set([str(x) for x in range(10)]))


class TestModelsStatic(TestCase):
    """Tests for static methods.
    """

    def test_clean_filesystem_location(self):
        FAKE_ABS_ROOT = '/root/of/all/evil'
        EXPECTED_CLEAN_URL = 'projects/blah'
        dirty_full_url = os.path.join(FAKE_ABS_ROOT, settings.MEDIA_ROOT,
                EXPECTED_CLEAN_URL)
        clean_location = clean_filesystem_location(dirty_full_url)
        self.assertEqual(EXPECTED_CLEAN_URL, clean_location)

    def test_auto_generate_short_name(self):
        LONG_NAME = 'E Coli K12 Substrain MG1655'
        EXPECTED_SHORT_NAME = 'e_coli_k12_s'
        self.assertEqual(EXPECTED_SHORT_NAME,
                auto_generate_short_name(LONG_NAME))
