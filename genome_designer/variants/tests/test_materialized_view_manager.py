"""
Tests for materialized_view_manager.py.
"""

from django.test import TestCase

from variants.materialized_view_manager import MATERIALIZED_TABLE_QUERYABLE_FIELDS_MAP


class TestMaterializedViewSchema(TestCase):

    def test_queryable_fields_map(self):
        """Make sure the schema builds without errors.
        """
        # Just check one of them and make sure it worked.
        position_schema = MATERIALIZED_TABLE_QUERYABLE_FIELDS_MAP['position']
        self.assertEquals(position_schema['type'], 'Integer')
