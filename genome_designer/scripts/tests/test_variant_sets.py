"""
Tests for adding and removing variants from variant_sets.
"""

import os

from django.contrib.auth.models import User
from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.bootstrap_data import create_fake_variants_and_variant_sets
from scripts.import_util import import_reference_genome_from_local_file
from scripts.variant_sets import add_or_remove_variants_from_set
import settings


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'
TEST_PROJECT_NAME = 'recoli'
REF_GENOME_1_LABEL = 'mg1655'
TEST_FASTA  = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')
SAMPLE_1_LABEL = 'sample1'
VARIANTSET_1_LABEL = 'Set A'
VARIANTSET_2_LABEL = 'Set B'


class TestAddAndRemoveVariantsFromSet(TestCase):
    """Tests for scripts.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        test_project = Project.objects.create(owner=user.get_profile(),
                title=TEST_PROJECT_NAME)
        self.ref_genome_1 = import_reference_genome_from_local_file(
            test_project, REF_GENOME_1_LABEL, TEST_FASTA, 'fasta')

        create_fake_variants_and_variant_sets(self.ref_genome_1)

        (self.sample_1, created) = ExperimentSample.objects.get_or_create(
                project=test_project,
                label=SAMPLE_1_LABEL)

        self.var_set1_uid = VariantSet.objects.get_or_create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_1_LABEL)[0].uid

        self.var_set2_uid = VariantSet.objects.get_or_create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_2_LABEL)[0].uid


    def test_add_variants_to_set(self):

        variant_uids = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__gt=25,
                chromosome='chrom').values_list('uid')

        response = add_or_remove_variants_from_set(
                variant_uids,
                'add',
                self.var_set1_uid)

        self.assertEqual(response['alert_type'], 'info', str(response))


    def test_remove_variant_from_set(self):

        variant_uids = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                chromosome='chrom',
                variantset__uid=self.var_set2_uid).values_list('uid')

        print self.var_set1_uid
        print variant_uids

        response = add_or_remove_variants_from_set(
                variant_uids,
                'remove',
                self.var_set2_uid)

        self.assertEqual(response['alert_type'], 'info', str(response))

