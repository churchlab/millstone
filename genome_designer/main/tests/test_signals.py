"""
Tests for signals.py.
"""

import os

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.model_utils import get_dataset_with_type

from scripts.import_util import import_reference_genome_from_local_file


TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'
TEST_PROJECT_NAME = 'testModels_project'
TEST_REF_GENOME_NAME = 'mg1655_partial'
TEST_REF_GENOME_PATH =  os.path.join(settings.PWD,
    'test_data/full_vcf_test_set/mg1655_tolC_through_zupT.gb')


class TestSignals(TestCase):

    def setUp(self):
        """Override.
        """
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)

        self.test_project = Project.objects.create(
            title=TEST_PROJECT_NAME,
            owner=user.get_profile())

        self.test_ref_genome = ReferenceGenome.objects.create(
            project=self.test_project,
            label='boom',
            num_chromosomes=2,
            num_bases=1000)

        self.test_ext_ref_genome = import_reference_genome_from_local_file(
            self.test_project,
            TEST_REF_GENOME_NAME,
            TEST_REF_GENOME_PATH,
            'genbank')

    def test_post_add_seq_to_ref_genome(self):
        """
        Ensure that everything gets converted after creating a new reference
        genome object, like snpeff, fasta, gff, etc.
        """

        # SNPEFF

        # check that the genbank file was symlinked
        gbk_path = os.path.join(
                self.test_ext_ref_genome.get_snpeff_directory_path(),
                'genes.gb')
        assert os.path.exists(gbk_path), 'snpeff gbk conversion failed: %s' % (
                gbk_path)

        # check that the db was made
        assert os.path.exists(os.path.join(
                self.test_ext_ref_genome.get_snpeff_directory_path(),
                'snpEffectPredictor.bin')), 'snpeff db was not made'

        # FASTA
        fasta = get_dataset_with_type(self.test_ext_ref_genome,
                type=Dataset.TYPE.REFERENCE_GENOME_FASTA)
        assert os.path.exists(fasta.get_absolute_location()), (
                'fasta conversion failed')

        # GFF
        gff = get_dataset_with_type(self.test_ext_ref_genome,
                type=Dataset.TYPE.REFERENCE_GENOME_GFF)
        assert os.path.exists(gff.get_absolute_location()), (
                'gff conversion failed')

    def test_post_variant_evidence_create(self):
        self.assertTrue(True)

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.test_ref_genome,
                chromosome='chrom',
                position=22,
                ref_value='A')

        va = VariantAlternate.objects.create(
                variant=variant,
                alt_value='G')

        sample_1 = ExperimentSample.objects.create(
                project=self.test_project,
                label='a sample')

        fake_source_dataset = Dataset.objects.create(
                type=Dataset.TYPE.VCF_FREEBAYES,
                label='fake',
                filesystem_location='')

        common_data_obj = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset=fake_source_dataset)

        # Test creating VE without data.
        ve_no_data = VariantEvidence.objects.create(
                experiment_sample=sample_1,
                variant_caller_common_data=common_data_obj)
        self.assertEqual(0, ve_no_data.variantalternate_set.all().count())

        # Test creating VE with ref gt.
        ve_ref_data = {
            'gt_bases': 'A/A'
        }
        ve_with_ref_data = VariantEvidence.objects.create(
                experiment_sample=sample_1,
                variant_caller_common_data=common_data_obj,
                data=ve_ref_data)
        self.assertEqual(0, ve_with_ref_data.variantalternate_set.all().count())

        # Test creating VE with alt gt.
        ve_alt_data = {
            'gt_bases': 'G/G'
        }
        ve_with_alt_data = VariantEvidence.objects.create(
                experiment_sample=sample_1,
                variant_caller_common_data=common_data_obj,
                data=ve_alt_data)
        self.assertEqual(1, ve_with_alt_data.variantalternate_set.all().count())
        self.assertEqual('G',
                ve_with_alt_data.variantalternate_set.all()[0].alt_value)
