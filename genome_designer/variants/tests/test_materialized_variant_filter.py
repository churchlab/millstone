"""
Tests for materialized_variant_filter.py.
"""

import os

from django.db import connection
from django.db import transaction
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
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.dynamic_snp_filter_key_map import update_filter_key_map
from settings import PWD as GD_ROOT
from variants.materialized_variant_filter import get_variants_that_pass_filter
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager


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

        self.materialized_view_manager = MeltedVariantMaterializedViewManager(
                self.ref_genome)


class TestVariantFilter(BaseTestVariantFilterTestCase):
    def test_filter__by_position(self):
        """Test filtering by position.
        """
        # Only variants that are part of a set will be returned in the melted
        # view. If a Variant is not associatd with a VariantSet, nor is it
        # associated with Sample, then the materialized view will not return
        # it.
        var_set = VariantSet.objects.create(
                reference_genome=self.ref_genome,
                label='set1')

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

            VariantToVariantSet.objects.create(variant=var, variant_set=var_set)

        # Test querying Variants with position > 5.
        result = get_variants_that_pass_filter('position > 5', self.ref_genome)
        variants_above_5 = result.variant_set
        self.assertEqual(4, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var['position'] > 5)

        # Test querying Variants with position >= 5.
        result = get_variants_that_pass_filter('position >= 5', self.ref_genome)
        variants_above_5 = result.variant_set
        self.assertEqual(5, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var['position'] >= 5)


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
