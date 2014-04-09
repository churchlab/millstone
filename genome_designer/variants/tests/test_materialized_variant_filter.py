"""
Tests for materialized_variant_filter.py.
"""

import os

from django.db import connection
from django.contrib.auth.models import User
from django.test import TestCase
from sympy.core.function import sympify

from main.models import AlignmentGroup
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
from main.testing_util import create_common_entities
from scripts.dynamic_snp_filter_key_map import MAP_KEY__ALTERNATE
from scripts.dynamic_snp_filter_key_map import update_filter_key_map
from settings import PWD as GD_ROOT
from variants.common import determine_visible_field_names
from variants.common import ParseError
from variants.materialized_variant_filter import get_variants_that_pass_filter
from variants.materialized_variant_filter import VariantFilterEvaluator
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__CHROMOSOME
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__POSITION


TEST_DIR = os.path.join(GD_ROOT, 'test_data', 'genbank_aligned')

TEST_ANNOTATED_VCF = os.path.join(TEST_DIR, 'bwa_align_annotated.vcf')


def run_query(filter_string, ref_genome):
    """Helper method that uses the VariantFilterEvaluator.
    """
    query_args = {
        'filter_string': filter_string,
        'visible_key_names': determine_visible_field_names(
            [], filter_string, ref_genome)
    }
    return get_variants_that_pass_filter(query_args, ref_genome)


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
        update_filter_key_map(self.ref_genome, TEST_ANNOTATED_VCF)

        self.vcf_dataset = Dataset.objects.create(
                label='test_data_set',
                type=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location=TEST_ANNOTATED_VCF)

        self.sample_obj_1 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample')

        self.sample_obj_2 = ExperimentSample.objects.create(
                project=self.project,
                label='fake sample 2')

        # Only variants that are part of a set will be returned in the melted
        # view. If a Variant is not associatd with a VariantSet, nor is it
        # associated with Sample, then the materialized view will not return
        # it.
        self.catchall_variant_set = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='catchall')

        self.materialized_view_manager = MeltedVariantMaterializedViewManager(
                self.ref_genome)


class TestVariantFilter(BaseTestVariantFilterTestCase):

    def test_filter__by_position(self):
        """Test filtering by position.
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

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        # Test querying Variants with position > 5.
        variants_above_5 = run_query('position > 5', self.ref_genome)
        self.assertEqual(4, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var[MELTED_SCHEMA_KEY__POSITION] > 5)

        # Test querying Variants with position >= 5.
        variants_above_5 = run_query('position >= 5', self.ref_genome)
        self.assertEqual(5, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var[MELTED_SCHEMA_KEY__POSITION] >= 5)


    def test_filter__by_chromosome(self):
        """Test filtering by chromosome.
        """
        CHROM_1 = 'chrom'
        for pos in range(6):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A')

            var.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=var,
                            alt_value='G'))

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        CHROM_2 = 'chrom2'
        for pos in range(9):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_2,
                position=pos,
                ref_value='A')

            var.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=var,
                            alt_value='G'))

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        def _assert_results(chromosome_value, num_results):
            variants = run_query(
                    'chromosome = %s'% chromosome_value, self.ref_genome)
            self.assertEqual(num_results, len(variants))
            for var in variants:
                self.assertEqual(chromosome_value, var[MELTED_SCHEMA_KEY__CHROMOSOME])
        _assert_results(CHROM_1, 6)
        _assert_results(CHROM_2, 9)


    def test_filter__by_position_and_chromosome(self):
        """Test filtering by position.
        """
        CHROM_1 = 'chrom'
        for pos in range(6):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A')

            var.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=var,
                            alt_value='G'))

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        CHROM_2 = 'chrom2'
        for pos in range(9):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_2,
                position=pos,
                ref_value='A')

            var.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=var,
                            alt_value='G'))

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)


        def _assert_results(query, num_results, chromosome_value,
                position_value, position_delim):
            variants = run_query(query, self.ref_genome)
            self.assertEqual(num_results, len(variants))
            for var in variants:
                self.assertEqual(chromosome_value, var[MELTED_SCHEMA_KEY__CHROMOSOME])
                self.assertTrue(eval(
                        str(var[MELTED_SCHEMA_KEY__POSITION]) + position_delim +
                                str(position_value)))

        QUERY_STRING = 'position > 4 & chromosome = chrom'
        _assert_results(QUERY_STRING, 1, CHROM_1, 4, '>')

        QUERY_STRING = 'position >= 5 & chromosome = chrom2'
        _assert_results(QUERY_STRING, 4, CHROM_2, 5, '>=')


    def test_filter__invalid_key(self):
        CHROM_1 = 'chrom'
        for pos in range(6):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A')

            var.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=var,
                            alt_value='G'))

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        QUERY_STRING = 'dinosaur > 4 & chromosome = chrom'
        with self.assertRaises(ParseError):
            run_query(QUERY_STRING, self.ref_genome)


    def test_filter__by_position_complex(self):
        """Test filtering by position.
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

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        # Test AND case.
        QUERY_STRING = 'position < 1 & position > 7'
        result = run_query(QUERY_STRING,
                self.ref_genome)
        self.assertEqual(0, len(result))

        # Test OR case.
        QUERY_STRING = 'position < 1 | position > 7'
        variant_list = run_query(QUERY_STRING,
                self.ref_genome)
        self.assertEqual(3, len(variant_list))
        for pos in [0, 8, 9]:
            found = False
            for var in variant_list:
                if var[MELTED_SCHEMA_KEY__POSITION] == pos:
                    found = True
                    break
            self.assertTrue(found, "Expected variant at pos %d not found" % pos)


    def test_filter__equality(self):
        """Test filtering with equality operators.
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

            VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        variants = run_query('position == 5',
                self.ref_genome)
        self.assertEqual(1, len(variants))
        self.assertEqual(5, variants.pop()[MELTED_SCHEMA_KEY__POSITION])

        variants = run_query('position = 5',
                self.ref_genome)
        self.assertEqual(1, len(variants))
        for var in variants:
            self.assertEqual(5, var[MELTED_SCHEMA_KEY__POSITION])

        variants = run_query('position != 5',
                self.ref_genome)
        self.assertEqual(9, len(variants))
        for var in variants:
            self.assertNotEqual(5, var[MELTED_SCHEMA_KEY__POSITION])


    def test_filter__no_sets(self):
        """Tests that we can filter on Variants that are not associated with
        sets. However, these Variants must be associated with ExperimentSamples.
        """
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

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

        common_data_obj = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset=self.vcf_dataset,
                alignment_group= alignment_group
        )

        raw_sample_data_dict = {
                'called': True,
                'gt_type': 2,
                'gt_bases': 'G/G'
        }

        sample_1_evidence = VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_1,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        raw_sample_data_dict = {
                'called': True,
                'gt_type': 2,
                'gt_bases': 'T/T'
        }
        sample_2_evidence = VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_2,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        variants = run_query('', self.ref_genome)

        # We expect 2 rows of results.
        self.assertEqual(2, len(variants))

        # Account for all results.
        sample_uid__alt_value_pairs = [
            (self.sample_obj_1.uid, 'G'),
            (self.sample_obj_2.uid, 'T'),
        ]


    def test_filter__key_missing_from_map(self):
        """Tests that filtering using an unrecognized key fails when using an
        improperly initialized key map.

        This test would catch something like this:
            https://github.com/churchlab/genome-designer-v2/issues/36
        """
        ref_genome_2 = ReferenceGenome.objects.create(project=self.project,
                label='refgenome2', num_chromosomes=1, num_bases=1000)

        # Initialize but don't update with source vcf, thus only global
        # keys are available.

        var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=ref_genome_2,
                chromosome='chrom',
                position=100,
                ref_value='A')

        var.variantalternate_set.add(
            VariantAlternate.objects.create(
                    variant=var,
                    alt_value='G'))

        VariantToVariantSet.objects.create(variant=var,
                    variant_set=self.catchall_variant_set)

        # This query runs without errors.
        QUERY_STRING = 'position < 1'
        variants = run_query(QUERY_STRING,
                ref_genome_2)
        self.assertEqual(0, len(variants))

        # This throws an error since INFO_XRM is not a recognized key.
        with self.assertRaises(ParseError):
            QUERY_STRING = 'position < 1 & INFO_XRM > 0'
            variants = run_query(QUERY_STRING,
                    ref_genome_2)


    def test_filter_by_json_field(self):
        """Filter using json fields.
        """
        # Make sure ref_genome_map has the field we are testing.
        alternate_key_map = self.ref_genome.variant_key_map[MAP_KEY__ALTERNATE]
        alternate_key_map['INFO_EFF_GENE'] = {
            u'num': -1,
            u'type': u'String'
        }
        self.ref_genome.save()

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        alt_data_dict = {
            'INFO_EFF_GENE': 'tolC'
        }
        VariantAlternate.objects.create(
                variant=variant,
                alt_value='T',
                data=alt_data_dict
        )

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

        common_data_obj = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset=self.vcf_dataset,
            alignment_group=alignment_group
        )

        raw_sample_data_dict = {
            'GT_BASES': 'T/T',
            'INFO_EFF_GENE': 'tolC'
        }
        VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_1,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        passing_variants = run_query('INFO_EFF_GENE = tolC',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))


    def test_case_insensitive(self):
        """Filter keys should not be case sensitive.

        Here we test a few. We might want more exhaustive tests.
        """
        # Fake the key map.
        alternate_key_map = self.ref_genome.variant_key_map[MAP_KEY__ALTERNATE]
        alternate_key_map['INFO_EFF_GENE'] = {
            u'num': -1,
            u'type': u'String'
        }
        self.ref_genome.save()

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=2,
                ref_value='A')

        alt_data_dict = {
            'INFO_EFF_GENE': 'tolC'
        }
        VariantAlternate.objects.create(
                variant=variant,
                alt_value='T',
                data=alt_data_dict
        )

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

        common_data_obj = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset=self.vcf_dataset,
            alignment_group=alignment_group)

        raw_sample_data_dict = {
            'GT_BASES': 'T/T',
            'GT_TYPE': 2,
            'INFO_EFF_GENE': 'tolC'
        }
        VariantEvidence.objects.create(
                experiment_sample=self.sample_obj_1,
                variant_caller_common_data=common_data_obj,
                data=raw_sample_data_dict)

        # Test empty query that should return everything.
        passing_variants = run_query('',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))

        # Test basic position query.
        passing_variants = run_query('position = 2',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))

        # Test both cases.
        passing_variants = run_query('POSITION = 2',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))

        # Test JSON field, uppercase.
        passing_variants = run_query('INFO_EFF_GENE = tolC',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))

        # Test JSON field, lowercase.
        passing_variants = run_query('info_eff_gene = tolC',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))

        # VariantEvidence data.
        passing_variants = run_query('gt_type = 2',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))
        passing_variants = run_query('GT_TYPE = 2',
                self.ref_genome)
        self.assertEqual(1, len(passing_variants))


class TestVariantFilterEvaluator(BaseTestVariantFilterTestCase):
    """Tests for the object that encapsulates evaluation of the filter string.
    """

    def test_variant_filter_constructor(self):
        """Tests the constructor.
        """
        query_args = {'filter_string': 'position > 5'}
        evaluator = VariantFilterEvaluator(query_args, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual('position > 5',
                evaluator.symbol_to_expression_map['A'])

        query_args = {'filter_string': 'position>5 & chromosome= chrom1'}
        evaluator = VariantFilterEvaluator(query_args, self.ref_genome)
        EXPECTED_SYMBOLIC_REP = sympify('A & B')
        self.assertEqual(EXPECTED_SYMBOLIC_REP, evaluator.sympy_representation)
        self.assertEqual('position>5',
                evaluator.symbol_to_expression_map['A'])
        self.assertEqual('chromosome= chrom1',
                evaluator.symbol_to_expression_map['B'])




class TestMinimal(BaseTestVariantFilterTestCase):
    """Minimal tests for materialized views.

    These were written in response to the initial set of tests not working, so
    I wanted to test a much simpler use of materialied views to rule out any
    issues potentially arising due to the Django testing flow's use of the
    database.
    """

    def test_minimal_materialized_view(self):
        """Creates a simple materialized view to check that these can be tested.
        """
        self.assertTrue(True)

        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A')

        # Sanity check.
        self.assertEqual(10,
                Variant.objects.filter(
                        reference_genome=self.ref_genome).count())

        cursor = connection.cursor()

        # Create a simple materialized view to check whether they even work.
        create_sql_statement = (
            'CREATE MATERIALIZED VIEW test_materialized_view AS '
                    'SELECT main_variant.id FROM main_variant '
        )
        cursor.execute(create_sql_statement)

        query_sql_statement = 'SELECT * FROM test_materialized_view'
        cursor.execute(query_sql_statement)

        self.assertEqual(10, len(cursor.fetchall()))


class TestDataConsistency(TestCase):
    """Tests data changing and the materialized view being correctly up to
    date.
    """

    def setUp(self):
        self.common_entities = create_common_entities()
        self.cursor = connection.cursor()

    def test_variant_set_membership_change(self):
        """Tests that the materialized view is refreshed when models change.
        """
        results = run_query('', self.common_entities['reference_genome'])
        self.assertEqual(0, len(results))

        # Add a Variant / VariantAlternate / VariantSet combo which we expect
        # to show up in the materialized view.
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.common_entities['reference_genome'],
                chromosome='chrom',
                position=2,
                ref_value='A')
        VariantAlternate.objects.create(
                variant=variant,
                alt_value='T')
        variant_set = VariantSet.objects.create(
                label='vs1',
                reference_genome=self.common_entities['reference_genome'])
        vtvs = VariantToVariantSet.objects.create(
                variant=variant,
                variant_set=variant_set)

        results = run_query('', self.common_entities['reference_genome'])
        self.assertEqual(1, len(results))

        # Delete the VariantSet association.
        vtvs.delete()

        results = run_query('', self.common_entities['reference_genome'])
        self.assertEqual(0, len(results))
