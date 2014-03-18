"""
Tests for melted_variant_schema.py.
"""

from django.test import TestCase

from variants.melted_variant_schema import MATERIALIZED_TABLE_QUERYABLE_FIELDS_MAP


class TestMaterializedViewSchema(TestCase):

    def test_queryable_fields_map(self):
        """Make sure the schema builds without errors.
        """
        # Just check one of them and make sure it worked.
        position_schema = MATERIALIZED_TABLE_QUERYABLE_FIELDS_MAP['position']
        self.assertEquals(position_schema['type'], 'Integer')
