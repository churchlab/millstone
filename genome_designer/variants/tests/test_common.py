"""
Tests for variants/common.py.
"""

import os

from django.contrib.auth.models import User
from django.test import TestCase

from main.models import Dataset
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantAlternate
from main.testing_util import create_common_entities_w_variants
from variants.dynamic_snp_filter_key_map import update_filter_key_map
from settings import PWD as GD_ROOT
from variants.common import determine_visible_field_names
from variants.common import extract_filter_keys
from variants.common import SymbolGenerator
from variants.common import update_parent_child_variant_fields


TEST_DIR = os.path.join(GD_ROOT, 'test_data', 'genbank_aligned')

TEST_ANNOTATED_VCF = os.path.join(TEST_DIR, 'bwa_align_annotated.vcf')


class TestCommon(TestCase):

    def setUp(self):
        user = User.objects.create_user('testuser', password='password',
                email='test@test.com')
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=self.project,
                label='refgenome')

        # Make sure the reference genome has the required vcf keys.
        update_filter_key_map(self.ref_genome, TEST_ANNOTATED_VCF)

        self.vcf_dataset = Dataset.objects.create(
                label='test_data_set',
                type=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location=TEST_ANNOTATED_VCF)


    def test_extract_filter_keys(self):
        """Tests extracting filter keys.
        """
        FILTER_EXPR = 'position > 5'
        EXPECTED_FILTER_KEY_SET = set(['POSITION'])
        self.assertEqual(EXPECTED_FILTER_KEY_SET,
                set(extract_filter_keys(FILTER_EXPR, self.ref_genome)))

        FILTER_EXPR = '(position < 5 & gt_type = 2) in ANY(1234, 4567)'
        EXPECTED_FILTER_KEY_SET = set(['POSITION', 'GT_TYPE'])
        self.assertEqual(EXPECTED_FILTER_KEY_SET,
                set(extract_filter_keys(FILTER_EXPR, self.ref_genome)))

    def test_determine_visible_field_names(self):
        EXPECTED_VISIBLE_KEYS = ['INFO_EFF_EFFECT']
        self.assertEqual(EXPECTED_VISIBLE_KEYS,
                determine_visible_field_names(EXPECTED_VISIBLE_KEYS, '',
                        self.ref_genome))


    def test_update_parent_child_variant_fields(self):
        self.common_entities = create_common_entities_w_variants()

        self.common_entities['samples'][0].add_child(
            self.common_entities['samples'][1])

        self.common_entities['samples'][0].add_child(
            self.common_entities['samples'][2])

        self.common_entities['samples'][2].add_child(
            self.common_entities['samples'][3])

        self.common_entities['samples'][4].add_child(
            self.common_entities['samples'][5])

        self.common_entities['samples'][5].add_child(
            self.common_entities['samples'][6])

        update_parent_child_variant_fields(
                self.common_entities['alignment_group'])

        v_1808 = Variant.objects.get(
                reference_genome=self.common_entities['reference_genome'],
                position=1808)
        self.assertEqual(set(v_1808.get_alternates()), set(['A']))

        vcc_1808 = v_1808.variantcallercommondata_set.get()
        ve_for_uid = lambda uid: (
            vcc_1808.variantevidence_set.get(experiment_sample__uid=uid))
        ve_sample_3 = ve_for_uid(u'9dd7a7a1')
        ve_sample_2 = ve_for_uid(u'9b19e708')

        self.assertEqual(ve_sample_3.data['GT_TYPE'],2)
        self.assertEqual(ve_sample_2.data['IN_CHILDREN'],1)


class TestSymbolGenerator(TestCase):
    """Tests the symbol generator used for symbolic manipulation.
    """

    def test_generator(self):
        symbol_maker = SymbolGenerator()
        self.assertEqual('A', symbol_maker.next())
        self.assertEqual('B', symbol_maker.next())
        self.assertEqual('C', symbol_maker.next())
