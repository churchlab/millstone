"""Tests for data_export_util.py.
"""

import StringIO

from django.test import TestCase
import vcf

from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities
from utils.data_export_util import export_variant_set_as_vcf


class TestExportVariantSetAsVcf(TestCase):

    def setUp(self):
        """Override.
        """
        self.common_entities = create_common_entities()

    def test_basic(self):
        variant_set = VariantSet.objects.create(
                reference_genome=self.common_entities['reference_genome'],
                label='vs1')

        for position in range(1, 11):
            var = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.common_entities['reference_genome'],
                    chromosome='chrom',
                    position=position,
                    ref_value='A')

            VariantAlternate.objects.create(
                    variant=var, alt_value='G')

            VariantToVariantSet.objects.create(
                    variant=var, variant_set=variant_set)

        output_fh = StringIO.StringIO()

        export_variant_set_as_vcf(variant_set, output_fh)

        output_fh.seek(0)

        reader = vcf.Reader(output_fh)
        row_count = 0
        for record in reader:
            self.assertEqual('G', record.ALT[0])
            row_count += 1

        self.assertEqual(10, row_count)
