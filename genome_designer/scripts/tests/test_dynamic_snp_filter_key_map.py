"""
Tests for dynamic_snp_filter_key_map.py.
"""

import os

from django.conf import settings
from django.test import TestCase

from main.testing_util import create_common_entities
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.dynamic_snp_filter_key_map import update_filter_key_map


TEST_VCF = os.path.join(settings.PWD, 'test_data',
        'fix_recoli_variants_snpeff_small_test.vcf')


class TestDynamicSnpFilterKeyMap(TestCase):

    def setUp(self):
        self.test_entities = create_common_entities()
        self.reference_genome = self.test_entities['reference_genome']

    def test_initialize(self):
        self.reference_genome.variant_key_map = initialize_filter_key_map()
        self.reference_genome.save()

    def test_update(self):
        self.reference_genome.variant_key_map = initialize_filter_key_map()
        self.reference_genome.save()
        update_filter_key_map(self.reference_genome, TEST_VCF)

    def test_update__no_initialize(self):
        update_filter_key_map(self.reference_genome, TEST_VCF)
