"""
Tests for adding and removing variants from variant_sets.
"""

import os

from django.conf import settings
from django.test import TestCase

from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantSet
from main.test_util import create_common_entities
from scripts.bootstrap_data import create_fake_variants_and_variant_sets
from variants.variant_sets import MODIFY_VARIANT_SET_MEMBERSHIP__ADD
from variants.variant_sets import MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE
from variants.variant_sets import update_variant_in_set_memberships


TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')
SAMPLE_1_LABEL = 'sample1'
VARIANTSET_1_LABEL = 'New Set A'
VARIANTSET_2_LABEL = 'New Set B'


class TestAddAndRemoveVariantsFromSet(TestCase):

    def setUp(self):
        common_entities = create_common_entities()
        project = common_entities['project']
        self.ref_genome_1 = common_entities['reference_genome']

        create_fake_variants_and_variant_sets(self.ref_genome_1)

        (self.sample_1, created) = ExperimentSample.objects.get_or_create(
                project=project,
                label=SAMPLE_1_LABEL)

        self.var_set1 = VariantSet.objects.create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_1_LABEL)

        self.var_set2 = VariantSet.objects.create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_2_LABEL)

    def test_add_and_remove(self):
        """Test add and remove without samples involved.

        TODO: Make generating Variants deterministic. Right now we use
        random positions.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        variant_obj_list = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__gt=25,
                chromosome='chrom')
        variant_uid_sample_uid_pair_list = [obj.uid
                for obj in variant_obj_list]

        ### Test add.
        response = update_variant_in_set_memberships(
                self.ref_genome_1,
                variant_uid_sample_uid_pair_list,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid)
        self.assertEqual(response['alert_type'], 'info', str(response))
        self.assertEqual(len(variant_obj_list),
                self.var_set1.variants.all().count())

        ### Test remove.
        num_variants_in_set1_before_remove = (
                self.var_set1.variants.all().count())
        self.assertTrue(self.var_set1.variants.all().count() > 0)
        variant_obj_list = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__gt=75,
                chromosome='chrom',
                variantset__uid=self.var_set1.uid)
        variant_uids_to_remove = [obj.uid for obj in variant_obj_list]
        self.assertTrue(len(variant_uids_to_remove) > 0)

        response = update_variant_in_set_memberships(
                self.ref_genome_1,
                variant_uids_to_remove,
                MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE,
                self.var_set1.uid)
        self.assertEqual(response['alert_type'], 'info', str(response))

        # Check that we've reduced the number of variants.
        self.assertEqual(
                num_variants_in_set1_before_remove -
                        len(variant_uids_to_remove),
                self.var_set1.variants.all().count())
