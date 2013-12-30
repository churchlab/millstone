"""
Tests for variant_filter.py.
"""

import os
import pickle

from django.contrib.auth.models import User
from django.test import TestCase
from sympy.core.function import sympify

from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.dynamic_snp_filter_key_map import update_filter_key_map
from settings import PWD as GD_ROOT
from variants.common import eval_variant_set_filter_expr
from variants.variant_filter import EXPRESSION_REGEX
from variants.variant_filter import SAMPLE_SCOPE_REGEX
from variants.variant_filter import SAMPLE_SCOPE_REGEX_NAMED
from variants.variant_filter import SET_REGEX
from variants.variant_filter import SET_REGEX_NAMED
from variants.variant_filter import get_variants_that_pass_filter
from variants.variant_filter import ParseError
from variants.variant_filter import VariantFilterEvaluator

TEST_DIR = os.path.join(GD_ROOT, 'test_data', 'genbank_aligned')

TEST_ANNOTATED_VCF = os.path.join(TEST_DIR, 'bwa_align_annotated.vcf')


class BaseTestVariantFilterTestCase(TestCase):
    """Base class for the test classes here.
    """

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

        self.sample_obj_1 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample',
                group='Plate 1',
                well='A01',
                num_reads=100,
        )

        self.sample_obj_2 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 2',
                group='Plate 1',
                well='A02',
                num_reads=100,
        )


class TestVariantFilter(BaseTestVariantFilterTestCase):


    def test_filter__with_non_SQL_key(self):
        """Tests handling of non-SQL key.

        NOTE: This is not actually testing filtering, just showing that we
        can handle the key.
        """
        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A')

            var.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=var,
                            alt_value='G'))

        QUERY_STRING = 'position < 1 & INFO_XRM > 0'
        variants = get_variants_that_pass_filter(QUERY_STRING,
                self.ref_genome).variant_set
        self.assertEqual(0, len(variants))

        QUERY_STRING = 'position < 1 | position >= 7 & INFO_XRM > 0'
        variants = get_variants_that_pass_filter(QUERY_STRING,
                self.ref_genome).variant_set
        self.assertEqual(1, len(variants))


    def test_filter__common_data(self):
        """Test with VariantCallerCommonData filtering.
        """
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        variant.variantalternate_set.add(
                VariantAlternate.objects.create(
                        variant=variant,
                        alt_value='G'))

        raw_data_dict = {
                'INFO_XRM': pickle.dumps(0.12),
        }
        VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset=self.vcf_dataset,
                data=raw_data_dict
        )

        QUERY_STRING = 'position < 5 & INFO_XRM > 0'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))

        QUERY_STRING = 'position < 5 & INFO_XRM > 1'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))


    def test_filter__variant_evidence(self):
        """Test with VariantEvidence objects.
        """
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        variant.variantalternate_set.add(
                VariantAlternate.objects.create(
                        variant=variant,
                        alt_value='G'))

        raw_common_data_dict = {
                'INFO_XRM': pickle.dumps(0.12),
        }
        common_data_obj = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset=self.vcf_dataset,
                data=raw_common_data_dict
        )

        sample_obj = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample',
                group='Plate 1',
                well='A01',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'gt_bases': pickle.dumps('G/G')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        sample_obj_2 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 2',
                group='Plate 1',
                well='A02',
                num_reads=100,
        )
        raw_sample_data_dict_2 = {
                'called': True,
                'gt_type': pickle.dumps(0),
                'gt_bases': pickle.dumps('A/A')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_2,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict_2)

        sample_obj_3 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 3',
                group='Plate 1',
                well='A03',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': False,
                'gt_type': pickle.dumps(None),
                'gt_bases': pickle.dumps(None)

        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_3,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)


        QUERY_STRING = 'position < 5 & gt_type = 2'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))

        QUERY_STRING = 'position < 5 & gt_type = 0'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))

        QUERY_STRING = 'position < 5 & gt_type = 1'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        QUERY_STRING = 'position > 5 & gt_type = 2'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        QUERY_STRING = 'position < 5 & (gt_type = 0 | gt_type = 2)'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))

        # Check that the number of passing samples is correct.
        passing_sample_ids = result.variant_id_to_metadata_dict[variant.id][
                'passing_sample_ids']
        self.assertEqual(2, len(passing_sample_ids))


    def test_filter__scoped(self):
        """Test with VariantEvidence objects.
        """
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        variant.variantalternate_set.add(
                VariantAlternate.objects.create(
                        variant=variant,
                        alt_value='G'))

        # Common data shared by all.
        raw_common_data_dict = {
                'INFO_XRM': pickle.dumps(0.12)
        }
        common_data_obj = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset=self.vcf_dataset,
                data=raw_common_data_dict
        )

        # Sample 1
        sample_obj = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample',
                group='Plate 1',
                well='A01',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'gt_bases': pickle.dumps('G/G')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        # Sample 2
        sample_obj_2 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 2',
                group='Plate 1',
                well='A02',
                num_reads=100,
        )
        raw_sample_data_dict_2 = {
                'called': True,
                'gt_type': pickle.dumps(0),
                'gt_bases': pickle.dumps('A/A')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_2,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict_2)

        # Sample 3
        sample_obj_3 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 3',
                group='Plate 1',
                well='A03',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': False,
                'gt_type': pickle.dumps(1),
                'gt_bases': pickle.dumps(None)
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_3,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        # Sample 4
        sample_obj_4 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 4',
                group='Plate 1',
                well='A04',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'gt_bases': pickle.dumps('G/G')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_4,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)


        QUERY_STRING = '(position < 5 & gt_type = 2) in ANY(%s, %s, %s)' % (
                sample_obj.uid, sample_obj_2.uid, sample_obj_3.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))

        QUERY_STRING = '(position < 5 & gt_type = 2) in ALL(%s, %s, %s)' % (
                sample_obj.uid, sample_obj_2.uid, sample_obj_3.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        QUERY_STRING = '(position < 5 & gt_type = 2) in ALL(%s)' % (
                sample_obj.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))

        QUERY_STRING = '(position < 5 & gt_type = 2) in ONLY(%s, %s, %s)' % (
                sample_obj.uid, sample_obj_2.uid, sample_obj_3.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        # This should fail because condition passes for sample_obj and
        # sample_obj_4.
        QUERY_STRING = '(position < 5 & gt_type = 2) in ONLY(%s)' % (
                sample_obj.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))


    def test_filter__common_data_per_alt(self):
        """Test filtering for common data of type '-1', 
        with different values for each alternate allele. 
        """

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        for alt_value, alt_eff in zip(['G','T'], ['NON_SYNONYMOUS_CODING',
                        'SYNONYMOUS_CODING']):

            alt = VariantAlternate.objects.create(
                        variant=variant,
                        alt_value=alt_value)
            alt.data = {'INFO_EFF_EFFECT': alt_eff}

            variant.variantalternate_set.add(alt)

        common_data_obj = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset=self.vcf_dataset,
                data=dict()
        )
        sample_obj = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample',
                group='Plate 1',
                well='A01',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'GT': pickle.dumps('1/1'),
                'gt_bases': pickle.dumps('G/G')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        sample_obj_2 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 2',
                group='Plate 1',
                well='A02',
                num_reads=100,
        )
        raw_sample_data_dict_2 = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'GT': pickle.dumps('1/2'),
                'gt_bases': pickle.dumps('G/T')
        }
        VariantEvidence.objects.create(
                uid=1000,
                experiment_sample=sample_obj_2,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict_2)

        sample_obj_3 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 3',
                group='Plate 1',
                well='A03',
                num_reads=100
        )
        raw_sample_data_dict = {
                'called': False,
                'gt_type': pickle.dumps(1),
                'GT': pickle.dumps('0/0'),
                'gt_bases': pickle.dumps('A/A')
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_3,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)


        # per_alt_dict, per_alt_types = get_per_alt_dict(
        #         'INFO_EFF_EFFECT', variant,
        #         VariantEvidence.objects.get(uid=1000),
        #         self.ref_genome.variant_key_map['snp_caller_common_data'])

        # assert('INFO_EFF_EFFECT' in per_alt_dict.keys())
        # assert('NON_SYNONYMOUS_CODING' in per_alt_dict['INFO_EFF_EFFECT'])

        QUERY_STRING = (
                '(position < 5 & INFO_EFF_EFFECT = NON_SYNONYMOUS_CODING)')
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set

        self.assertEqual(1, len(variants))  

        passing_sample_ids = result.variant_id_to_metadata_dict[variant.id][
                'passing_sample_ids']

        self.assertEqual(2, len(passing_sample_ids))

    def test_filter__sets(self):
        """Test filtering relative to sets.
        """
        # We create 3 VariantSets.
        set_1 = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set1')
        set_2 = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set2')
        set_3 = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set3')

        # Add some variants to each of the above sets.
        def _create_variants_in_set(pos_range, var_set_list):
            for pos in pos_range:
                var = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome,
                    chromosome='chrom',
                    position=pos,
                    ref_value='A')

                var.variantalternate_set.add(
                        VariantAlternate.objects.create(
                                variant=var,
                                alt_value='G'))

                for var_set in var_set_list:
                    if var_set is not None:
                        VariantToVariantSet.objects.create(variant=var,
                                variant_set=var_set)
        _create_variants_in_set(range(3), [set_1]) # 3
        _create_variants_in_set(range(3, 6), [set_2]) # 3
        _create_variants_in_set(range(6, 9), [set_3]) # 3
        _create_variants_in_set(range(9, 13), [None]) # 4
        _create_variants_in_set(range(13, 19), [set_1, set_2]) # 6

        # We have 19 total variants.
        self.assertEqual(19, len(Variant.objects.all()))

        # 3 + 6 = 9 in set_1
        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(9, len(variants))

        # 3 + 6 + 3 = 12 between set_1 and set_2
        QUERY_STRING = 'IN_SET(%s) | IN_SET(%s)' % (set_1.uid, set_2.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(12, len(variants))

        # All but 4, so 15, in any set.
        QUERY_STRING = 'IN_SET(%s) | IN_SET(%s) | IN_SET(%s)' % (
                set_1.uid, set_2.uid, set_3.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(15, len(variants))

        # All but 9, or 10 not in set_1.
        QUERY_STRING = 'NOT_IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(10, len(variants))

        QUERY_STRING = 'IN_SET(%s) & IN_SET(%s)' % (set_1.uid, set_3.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        QUERY_STRING = 'IN_SET(%s) & IN_SET(%s)' % (set_1.uid, set_2.uid)

        # Sanity check.
        var_in_set_1 = set(Variant.objects.filter(
                varianttovariantset__variant_set__uid=set_1.uid))
        var_in_set_2 = set(Variant.objects.filter(
                varianttovariantset__variant_set__uid=set_2.uid))
        self.assertEqual(6, len(var_in_set_1 & var_in_set_2))

        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(6, len(variants))


    def test_filter__sets_per_sample(self):
        """Test filtering over sets when associated with samples.
        """
        set_1 = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set1')

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=1,
                ref_value='A')

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        common_data_obj = VariantCallerCommonData.objects.create(
                variant=var,
                source_dataset=self.vcf_dataset)

        sample_1_evidence = VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_1,
                variant_caller_common_data=common_data_obj)

        vtvs = VariantToVariantSet.objects.create(
                variant=var, variant_set=set_1)

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set(),
                metadata[passing_variant.id]['passing_sample_ids'])

        # Add a sample association.
        vtvs.sample_variant_set_association.add(self.sample_obj_1)
        vtvs.sample_variant_set_association.add(self.sample_obj_2)

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([self.sample_obj_1.id]),
                metadata[passing_variant.id]['passing_sample_ids'])


    def test_filter__scoped__per_alt(self):
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        alt_G = VariantAlternate.objects.create(
                variant=variant,
                alt_value='G')
        variant.variantalternate_set.add(alt_G)

        alt_T = VariantAlternate.objects.create(
                variant=variant,
                alt_value='T')
        variant.variantalternate_set.add(alt_T)

        # Asserts before adding samples.
        query_and_num_expected_pairs = [
                ('alt_value = A', 0),
                ('alt_value = T', 1),
                ('alt_value = G', 1)
        ]
        for query_string, num_expected in query_and_num_expected_pairs:
            result = get_variants_that_pass_filter(query_string, self.ref_genome)
            variants = result.variant_set
            self.assertEqual(num_expected, len(variants),
                    "Expected: %d, Actual: %d, Query: %s" % (
                            num_expected, len(variants), query_string))

        # Create a common data object and samples with different alt_values.

        common_data_obj = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset=self.vcf_dataset,
        )

        sample_obj_1 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample',
                group='Plate 1',
                well='A01',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'gt_bases': pickle.dumps('G/G')
        }
        sample_1_evidence = VariantEvidence.objects.create(
                experiment_sample=sample_obj_1,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        sample_obj_2 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 2',
                group='Plate 1',
                well='A02',
                num_reads=100,
        )
        raw_sample_data_dict = {
                'called': True,
                'gt_type': pickle.dumps(2),
                'gt_bases': pickle.dumps('T/T')
        }
        sample_2_evidence = VariantEvidence.objects.create(
                experiment_sample=sample_obj_2,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        QUERY_STRING = 'alt_value = T'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([sample_obj_2.id]),
                metadata[passing_variant.id]['passing_sample_ids'])

        QUERY_STRING = 'alt_value != T'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([sample_obj_1.id]),
                metadata[passing_variant.id]['passing_sample_ids'])


class TestVariantFilterEvaluator(BaseTestVariantFilterTestCase):
    """Tests for the object that encapsulates evaluation of the filter string.
    """

    def test_variant_filter_constructor__scoped(self):
        """Tests the constructor.
        """
        raw_filter_string = '(position > 5) in ALL(s1)'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual(raw_filter_string,
                evaluator.symbol_to_expression_map['A'])

        raw_filter_string = '((position > 5) in ALL(s1))'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('(A)')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual(raw_filter_string,
                evaluator.symbol_to_expression_map['A'])

        raw_filter_string = '((position > 5) in ALL(s1)) | position < 2'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('(A) | B')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual('((position > 5) in ALL(s1))',
                evaluator.symbol_to_expression_map['A'])
        self.assertEqual('position < 2',
                evaluator.symbol_to_expression_map['B'])


    def test_variant_filter_constructor__sets(self):
        """Test parsing out the set annotation.
        """
        raw_filter_string = 'IN_SET(1234)'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation,
                "Expected: %s. Actual: %s" % (EXPECTED_SYMBOLIC_REP,
                        evaluator.sympy_representation))
        self.assertEqual('IN_SET(1234)',
                evaluator.symbol_to_expression_map['A'])


    def test_variant_filter_constructor__sets_combo(self):
        raw_filter_string = 'IN_SET(1234) | IN_SET(5678)'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A | B')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation,
                "Expected: %s. Actual: %s" % (EXPECTED_SYMBOLIC_REP,
                        evaluator.sympy_representation))
        self.assertEqual('IN_SET(1234)',
                evaluator.symbol_to_expression_map['A'])
        self.assertEqual('IN_SET(5678)',
                evaluator.symbol_to_expression_map['B'])

        raw_filter_string = 'IN_SET(1234) & IN_SET(5678)'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A & B')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation,
                "Expected: %s. Actual: %s" % (EXPECTED_SYMBOLIC_REP,
                        evaluator.sympy_representation))


class TestExpressionRegex(TestCase):
    """Tests the regular expression that recognizes key op value conditions.
    """

    def test_expression(self):
        """Basic expression tests.
        """
        # Test positive matches.
        self.assertTrue(EXPRESSION_REGEX.match('key=value'))
        self.assertTrue(EXPRESSION_REGEX.match('key==value'))
        self.assertTrue(EXPRESSION_REGEX.match('key!=value'))

        # Test positive, with spaces.
        self.assertTrue(EXPRESSION_REGEX.match('key = value'))

        # Test negative matches.
        self.assertFalse(EXPRESSION_REGEX.match('key = = value'))
        self.assertFalse(EXPRESSION_REGEX.match('keyvalue'))

        # Test complex matches.
        QUERY = '(gt_type = 2 | gt_type = 1 | position > 100) in ALL(fa3a8451)'
        self.assertTrue(SAMPLE_SCOPE_REGEX.match(QUERY))


    def test_sample_scope_regex(self):
        """Basic expression tests.
        """
        # Test positive matches.
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in ALL(bar)'))
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in ANY(bar)'))
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in ONLY(bar)'))
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in all(bar)'))
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in any(bar)'))
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in only(bar)'))
        self.assertTrue(SAMPLE_SCOPE_REGEX.match('(foo) in ALL(bar,two,three)'))

        # Test negative matches.
        self.assertFalse(SAMPLE_SCOPE_REGEX.match('(foo) in ALIEN(bar)'))
        self.assertFalse(SAMPLE_SCOPE_REGEX.match('(foo) in A(bar)'))
        self.assertFalse(SAMPLE_SCOPE_REGEX.match('foo in A(bar)'))

        # Test complex matches.
        QUERY = '(gt_type = 2 | gt_type = 1 | position > 100) in ALL(fa3a8451)'
        self.assertTrue(SAMPLE_SCOPE_REGEX.match(QUERY))


    def test_sample_scope_regex__named(self):
        QUERY_STRING = '(position > 5) in ALL(s1,s2)'
        match = SAMPLE_SCOPE_REGEX_NAMED.match(QUERY_STRING)
        self.assertEqual('position > 5', match.group('condition'))
        self.assertEqual('ALL', match.group('scope_type'))
        self.assertEqual('s1,s2', match.group('samples'))

        QUERY_STRING = '(position > 5 | gt_type = 2) in ONLY(s1)'
        match = SAMPLE_SCOPE_REGEX_NAMED.match(QUERY_STRING)
        self.assertEqual('position > 5 | gt_type = 2', match.group('condition'))
        self.assertEqual('ONLY', match.group('scope_type'))
        self.assertEqual('s1', match.group('samples'))


    def test_set_regex(self):
        QUERY_STRING = 'IN_SET(1234)'
        self.assertTrue(SET_REGEX.match(QUERY_STRING))

        QUERY_STRING = 'IN_SET(5678abda)'
        self.assertTrue(SET_REGEX.match(QUERY_STRING))

        QUERY_STRING = 'NOT_IN_SET(1234)'
        self.assertTrue(SET_REGEX.match(QUERY_STRING))

        # Match only the first part.
        QUERY_STRING = 'IN_SET(1234) | IN_SET(5678)'
        self.assertEqual('IN_SET(1234)', SET_REGEX.match(QUERY_STRING).group())


    def test_set_regex__named(self):
        QUERY_STRING = 'IN_SET(1234)'
        match = SET_REGEX_NAMED.match(QUERY_STRING)
        self.assertEqual('1234', match.group('sets'))
        self.assertFalse(match.group('maybe_not'))

        QUERY_STRING = 'NOT_IN_SET(5678)'
        match = SET_REGEX_NAMED.match(QUERY_STRING)
        self.assertEqual('5678', match.group('sets'))
        self.assertTrue(match.group('maybe_not'))


class TestFilterVariantSets(BaseTestVariantFilterTestCase):
    """Tests for filtering relative to VariantSets.
    """

    def test_filter_sets(self):
        self.assertTrue(True)

        # 2 different VariantSets.
        set_1 = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set1')

        set_2 = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set2')

        # 1 Variant.
        var = Variant.objects.create(
                 type=Variant.TYPE.TRANSITION,
                 reference_genome=self.ref_genome,
                 chromosome='chrom',
                 position=1,
                 ref_value='A')

        # Not added to set yet, so no results expected.
        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        # Now we add the Variant to the set, but don't associate it with a
        # sample yet.
        common_data_obj = VariantCallerCommonData.objects.create(
                variant=var,
                source_dataset=self.vcf_dataset)

        sample_1_evidence = VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_1,
                variant_caller_common_data=common_data_obj)

        vtvs = VariantToVariantSet.objects.create(
                variant=var, variant_set=set_1)

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set(),
                metadata[passing_variant.id]['passing_sample_ids'])

        # Add a sample association.
        vtvs.sample_variant_set_association.add(self.sample_obj_1)

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([self.sample_obj_1.id]),
                metadata[passing_variant.id]['passing_sample_ids'])

        # Add another sample association.
        sample_2_evidence = VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_2,
                variant_caller_common_data=common_data_obj)

        vtvs.sample_variant_set_association.add(self.sample_obj_2)

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([self.sample_obj_1.id, self.sample_obj_2.id]),
                metadata[passing_variant.id]['passing_sample_ids'])

        # Test the negative case.
        # Both Samples have the Variant in Set1 so there should be no results.
        QUERY_STRING = 'NOT_IN_SET(%s)' % (set_1.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(0, len(variants))

        # Now we test the NOT IN case when there is no association, use set_2.
        QUERY_STRING = 'NOT_IN_SET(%s)' % (set_2.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([self.sample_obj_1.id, self.sample_obj_2.id]),
                metadata[passing_variant.id]['passing_sample_ids'])

        # Now create an association with just one sample and test it out.
        vtvs_set_2 = VariantToVariantSet.objects.create(
                variant=var, variant_set=set_2)
        vtvs_set_2.sample_variant_set_association.add(self.sample_obj_1)

        QUERY_STRING = 'NOT_IN_SET(%s)' % (set_2.uid)
        result = eval_variant_set_filter_expr(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))
        passing_variant = list(variants)[0]
        metadata = result.variant_id_to_metadata_dict
        self.assertEqual(set([self.sample_obj_2.id]),
                metadata[passing_variant.id]['passing_sample_ids'])
