"""
Tests for materialized_view_manager.py.
"""

from django.db import connection
from django.test import TestCase

from main.models import Chromosome
from main.models import Dataset
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VS_LABEL
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VS_UID
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
                chromosome=Chromosome.objects.get(reference_genome=self.common_entities['reference_genome']),
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


    def test_with_samples__no_associations(self):
        """Test example with sample associations.
        """
        mvm = MeltedVariantMaterializedViewManager(
                self.common_entities['reference_genome'])
        # Both variant sets should be in same row, aggregated.
        mvm.create()
        self.cursor.execute('SELECT * FROM %s' % mvm.get_table_name())
        self.assertEqual(0, len(self.cursor.fetchall()))

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.common_entities['reference_genome'],
                chromosome=Chromosome.objects.get(reference_genome=self.common_entities['reference_genome']),
                position=2,
                ref_value='A'
        )
        VariantAlternate.objects.create(
                variant=variant,
                alt_value='T',
        )

        vcf_source_dataset = Dataset.objects.create(
            type=Dataset.TYPE.VCF_FREEBAYES,
            label='fake_source_dataset')

        common_data_obj = VariantCallerCommonData.objects.create(
                alignment_group=self.common_entities['alignment_group_1'],
                variant=variant,
                source_dataset=vcf_source_dataset)

        VariantEvidence.objects.create(
                experiment_sample=self.common_entities['sample_1'],
                variant_caller_common_data=common_data_obj,
        )

        VariantEvidence.objects.create(
                experiment_sample=self.common_entities['sample_2'],
                variant_caller_common_data=common_data_obj,
        )

        # Associate a VariantSet, but no sample association.
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
        results = [dict(zip([col[0].upper() for col in self.cursor.description], row))
                for row in self.cursor.fetchall()]

        # Expect 3 rows: 2 sample-associated, one catch-all.
        self.assertEqual(3, len(results))

        # Only 1, the catch-all row, should have VariantSet data.
        observed_rows_with_variant_set_data = 0
        for data_row in results:
            self.assertEqual(1, len(data_row[MELTED_SCHEMA_KEY__VS_UID]))
            if data_row[MELTED_SCHEMA_KEY__VS_UID][0] is not None:
                observed_rows_with_variant_set_data += 1
        self.assertEqual(1, observed_rows_with_variant_set_data)
