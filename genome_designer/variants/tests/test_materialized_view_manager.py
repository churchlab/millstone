"""
Tests for materialized_view_manager.py.
"""

from django.db import connection
from django.test import TestCase

from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager


class TestMaterializedViewManager(TestCase):
    """Tests for the variants materialized view.
    """

    def setUp(self):
        self.common_entities = create_common_entities()

        self.cursor = connection.cursor()

    def test_multiple_variant_sets(self):
        mvm = MeltedVariantMaterializedViewManager(
                self.common_entities['reference_genome'])

        # Expect 0 before adding any entities.
        mvm.create()
        self.cursor.execute('SELECT * FROM %s' % mvm.get_table_name())
        self.assertEqual(0, len(self.cursor.fetchall()))

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.common_entities['reference_genome'],
                chromosome='chrom',
                position=2,
                ref_value='A'
        )
        VariantAlternate.objects.create(
                variant=variant,
                alt_value='T',
        )

        # First variant set.
        variant_set = VariantSet.objects.create(
                label='vs1',
                reference_genome=self.common_entities['reference_genome']
        )
        VariantToVariantSet.objects.create(
                variant=variant,
                variant_set=variant_set
        )

        mvm.create()
        self.cursor.execute('SELECT * FROM %s' % mvm.get_table_name())
        self.assertEqual(1, len(self.cursor.fetchall()))

        # Second variant set.
        variant_set_2 = VariantSet.objects.create(
                label='vs2',
                reference_genome=self.common_entities['reference_genome']
        )
        VariantToVariantSet.objects.create(
                variant=variant,
                variant_set=variant_set_2
        )

        # Both variant sets should be in same row, aggregated.
        mvm.create()
        self.cursor.execute('SELECT * FROM %s' % mvm.get_table_name())
        results = self.cursor.fetchall()
        self.assertEqual(1, len(results))
