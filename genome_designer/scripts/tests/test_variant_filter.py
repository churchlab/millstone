"""
Tests for variant_filter.py.
"""

import pickle

from django.contrib.auth.models import User
from django.test import TestCase
from sympy.core.function import sympify

from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.variant_filter import EXPRESSION_REGEX
from scripts.variant_filter import SAMPLE_SCOPE_REGEX
from scripts.variant_filter import SAMPLE_SCOPE_REGEX_NAMED
from scripts.variant_filter import SET_REGEX
from scripts.variant_filter import SET_REGEX_NAMED
from scripts.variant_filter import get_variants_that_pass_filter
from scripts.variant_filter import symbol_generator
from scripts.variant_filter import ParseError
from scripts.variant_filter import VariantFilterEvaluator


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
        self.vcf_dataset = Dataset.objects.create(
                label='test_data_set',
                type=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location='fake location')



class TestVariantFilter(BaseTestVariantFilterTestCase):

    def test_filter__by_position(self):
        """Test filtering by position.
        """
        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A',
                alt_value='G')

        # Test querying Variants with position > 5.
        result = get_variants_that_pass_filter('position > 5', self.ref_genome)
        variants_above_5 = result.variant_set
        self.assertEqual(4, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var.position > 5)

        # Test querying Variants with position >= 5.
        result = get_variants_that_pass_filter('position >= 5', self.ref_genome)
        variants_above_5 = result.variant_set
        self.assertEqual(5, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var.position >= 5)


    def test_filter__by_chromosome(self):
        """Test filtering by position.
        """
        CHROM_1 = 'chrom'
        for pos in range(6):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A',
                alt_value='G')

        CHROM_2 = 'chrom2'
        for pos in range(9):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_2,
                position=pos,
                ref_value='A',
                alt_value='G')

        result = get_variants_that_pass_filter('chromosome = chrom',
                self.ref_genome)
        self.assertEqual(6, len(result.variant_set))

        result = get_variants_that_pass_filter('chromosome = chrom2',
                self.ref_genome)
        self.assertEqual(9, len(result.variant_set))


    def test_filter__by_position_and_chromosome(self):
        """Test filtering by position.
        """
        CHROM_1 = 'chrom'
        for pos in range(6):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A',
                alt_value='G')

        CHROM_2 = 'chrom2'
        for pos in range(9):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_2,
                position=pos,
                ref_value='A',
                alt_value='G')

        QUERY_STRING = 'position > 4 & chromosome = chrom'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        self.assertEqual(1, len(result.variant_set))

        QUERY_STRING = 'position >= 5 & chromosome = chrom2'
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        self.assertEqual(4, len(result.variant_set))


    def test_filter__invalid_key(self):
        CHROM_1 = 'chrom'
        for pos in range(6):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A',
                alt_value='G')

        QUERY_STRING = 'dinosaur > 4 & chromosome = chrom'
        with self.assertRaises(ParseError):
            get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)


    def test_filter__by_position_complex(self):
        """Test filtering by position.
        """
        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A',
                alt_value='G')

        # Test AND case.
        QUERY_STRING = 'position < 1 & position > 7'
        result = get_variants_that_pass_filter(QUERY_STRING,
                self.ref_genome)
        self.assertEqual(0, len(result.variant_set))

        # Test OR case.
        QUERY_STRING = 'position < 1 | position > 7'
        result = get_variants_that_pass_filter(QUERY_STRING,
                self.ref_genome)
        variant_list = result.variant_set
        self.assertEqual(3, len(variant_list))
        for pos in [0, 8, 9]:
            found = False
            for var in variant_list:
                if var.position == pos:
                    found = True
                    break
            self.assertTrue(found, "Expected variant at pos %d not found" % pos)


    def test_filter__with_non_SQL_key(self):
        """Tests handling of non-SQL key.

        NOTE: This is not actually testing filtering, just showing that we
        can handle the key.
        """
        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A',
                alt_value='G')

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
                ref_value='A',
                alt_value='G')
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
                ref_value='A',
                alt_value='G')

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
                'gt_type': 1,
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


    def test_filter__equality(self):
        """Test filtering with equality operators.
        """
        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A',
                alt_value='G')

        variants = get_variants_that_pass_filter('position == 5',
                self.ref_genome).variant_set
        self.assertEqual(1, len(variants))
        self.assertEqual(5, variants.pop().position)

        variants = get_variants_that_pass_filter('position = 5',
                self.ref_genome).variant_set
        self.assertEqual(1, len(variants))
        self.assertEqual(5, variants.pop().position)

        variants = get_variants_that_pass_filter('position != 5',
                self.ref_genome).variant_set
        self.assertEqual(9, len(variants))
        for var in variants:
            self.assertTrue(var.position != 5)


    def test_filter__scoped(self):
        """Test with VariantEvidence objects.
        """
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A',
                alt_value='G')

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
                'gt_type': pickle.dumps(1),
        }
        VariantEvidence.objects.create(
                experiment_sample=sample_obj_3,
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

        QUERY_STRING = '(position < 5 & gt_type = 2) in ONLY(%s)' % (
                sample_obj.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(1, len(variants))


    def test_filter__sets(self):
        """Test filtering relative to sets.
        """
        # Some variant sets for testing.
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
                    ref_value='A',
                    alt_value='G')
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

        QUERY_STRING = 'IN_SET(%s)' % (set_1.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(9, len(variants))

        QUERY_STRING = 'IN_SET(%s) | IN_SET(%s)' % (set_1.uid, set_2.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(12, len(variants))

        QUERY_STRING = 'IN_SET(%s) | IN_SET(%s) | IN_SET(%s)' % (
                set_1.uid, set_2.uid, set_3.uid)
        result = get_variants_that_pass_filter(QUERY_STRING, self.ref_genome)
        variants = result.variant_set
        self.assertEqual(15, len(variants))

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


class TestVariantFilterEvaluator(BaseTestVariantFilterTestCase):
    """Tests for the object that encapsulates evaluation of the filter string.
    """

    def test_variant_filter_constructor(self):
        """Tests the constructor.
        """
        raw_filter_string = 'position > 5'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual('position > 5',
                evaluator.symbol_to_expression_map['A'])

        raw_filter_string = 'position>5 & chromosome= chrom1'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A & B')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual('position>5',
                evaluator.symbol_to_expression_map['A'])
        self.assertEqual('chromosome= chrom1',
                evaluator.symbol_to_expression_map['B'])


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


class TestSymbolGenerator(TestCase):
    """Tests the symbol generator used for symbolic manipulation.
    """

    def test_generator(self):
        symbol_maker = symbol_generator()
        self.assertEqual('A', symbol_maker.next())
        self.assertEqual('B', symbol_maker.next())
        self.assertEqual('C', symbol_maker.next())


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
