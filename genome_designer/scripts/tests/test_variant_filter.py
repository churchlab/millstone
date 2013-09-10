"""
Tests for variant_filter.py.
"""

from django.contrib.auth.models import User
from django.test import TestCase
from sympy.core.function import sympify

from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
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
        project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=project,
                label='refgenome', num_chromosomes=1, num_bases=1000)


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
        variants_above_5 = get_variants_that_pass_filter('position > 5',
                self.ref_genome)
        self.assertEqual(4, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var.position > 5)

        # Test querying Variants with position >= 5.
        variants_above_5 = get_variants_that_pass_filter('position >= 5',
                self.ref_genome)
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

        self.assertEqual(6, len(get_variants_that_pass_filter(
                'chromosome = chrom', self.ref_genome)))

        self.assertEqual(9, len(get_variants_that_pass_filter(
                'chromosome = chrom2', self.ref_genome)))


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
        self.assertEqual(1, len(get_variants_that_pass_filter(
                QUERY_STRING, self.ref_genome)))

        QUERY_STRING = 'position >= 5 & chromosome = chrom2'
        self.assertEqual(4, len(get_variants_that_pass_filter(
                QUERY_STRING, self.ref_genome)))


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
        self.assertEqual('position>5',
                evaluator.symbol_to_expression_map['A'])

        raw_filter_string = 'position>5 & chromosome= chrom1'
        evaluator = VariantFilterEvaluator(raw_filter_string, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A & B')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual('position>5',
                evaluator.symbol_to_expression_map['A'])
        self.assertEqual('chromosome=chrom1',
                evaluator.symbol_to_expression_map['B'])


class TestSymbolGenerator(TestCase):
    """Tests the symbol generator used for symbolic manipulation.
    """

    def test_generator(self):
        symbol_maker = symbol_generator()
        self.assertEqual('A', symbol_maker.next())
        self.assertEqual('B', symbol_maker.next())
        self.assertEqual('C', symbol_maker.next())
