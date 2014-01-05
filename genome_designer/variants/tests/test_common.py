"""
Tests for variants/common.py.
"""

import os

from django.contrib.auth.models import User
from django.test import TestCase

from main.models import Dataset
from main.models import Project
from main.models import ReferenceGenome
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.dynamic_snp_filter_key_map import update_filter_key_map
from settings import PWD as GD_ROOT
from variants.common import extract_filter_keys
from variants.common import SymbolGenerator

TEST_DIR = os.path.join(GD_ROOT, 'test_data', 'genbank_aligned')

TEST_ANNOTATED_VCF = os.path.join(TEST_DIR, 'bwa_align_annotated.vcf')


class TestCommon(TestCase):

    def setUp(self):
        user = User.objects.create_user('testuser', password='password',
                email='test@test.com')
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=self.project,
                label='refgenome', num_chromosomes=1, num_bases=1000)

        # Make sure the reference genome has the required vcf keys.
        initialize_filter_key_map(self.ref_genome)
        update_filter_key_map(self.ref_genome, TEST_ANNOTATED_VCF)

        self.vcf_dataset = Dataset.objects.create(
                label='test_data_set',
                type=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location=TEST_ANNOTATED_VCF)


    def test_extract_filter_keys(self):
        """Tests extracting filter keys.
        """
        FILTER_EXPR = 'position > 5'
        EXPECTED_FILTER_KEY_SET = set(['position'])
        self.assertEqual(EXPECTED_FILTER_KEY_SET,
                set(extract_filter_keys(FILTER_EXPR, self.ref_genome)))

        FILTER_EXPR = '(position < 5 & gt_type = 2) in ANY(1234, 4567)'
        EXPECTED_FILTER_KEY_SET = set(['position', 'gt_type'])
        self.assertEqual(EXPECTED_FILTER_KEY_SET,
                set(extract_filter_keys(FILTER_EXPR, self.ref_genome)))


class TestSymbolGenerator(TestCase):
    """Tests the symbol generator used for symbolic manipulation.
    """

    def test_generator(self):
        symbol_maker = SymbolGenerator()
        self.assertEqual('A', symbol_maker.next())
        self.assertEqual('B', symbol_maker.next())
        self.assertEqual('C', symbol_maker.next())
