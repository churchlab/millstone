"""
Tests for model_views.py.
"""

from django.test import TestCase

from main.model_views import _modify_obj_list_for_variant_set_display
from main.model_views import ALL_VS_LABEL_KEY
from main.model_views import ALL_VS_UID_KEY
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VS_LABEL
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VS_UID


class TestModelViews(TestCase):

    def test_modify_obj_list_for_variant_set_display(self):
        """Test a catch-all row and a sample-associated row.
        """
        initial_obj_list = [
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_ID': None,
                    'EXPERIMENT_SAMPLE_LABEL': None,
                    'EXPERIMENT_SAMPLE_UID': None,
                    'VARIANT_SET_LABEL': ['green'],
                    'VARIANT_SET_UID': ['88349a95'],
                    'ALT': u'T',
                    'REF': u'A',
                },
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_ID': 'es1_id',
                    'EXPERIMENT_SAMPLE_LABEL': 'es1',
                    'EXPERIMENT_SAMPLE_UID': 'es1_uid',
                    'VARIANT_SET_LABEL': [],
                    'VARIANT_SET_UID': [],
                    'ALT': u'T',
                    'REF': u'A',
                }
        ]

        modified_obj_list = _modify_obj_list_for_variant_set_display(
                initial_obj_list)

        # Rows preserved.
        self.assertEqual(1, len(modified_obj_list))

        # Require keys preserved.
        for obj in modified_obj_list:
            self.assertEqual(1, len(obj['all_variant_set_label']))
            self.assertEqual(1, len(obj['all_variant_set_uid']))

    def test_modify_obj_list_for_variant_set_display__different_alts(self):
        """Test that we preserve the rows for different alts.
        """
        initial_obj_list = [
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_ID': None,
                    'EXPERIMENT_SAMPLE_LABEL': None,
                    'EXPERIMENT_SAMPLE_UID': None,
                    'VARIANT_SET_LABEL': ['green'],
                    'VARIANT_SET_UID': ['88349a95'],
                    'ALT': u'T',
                    'REF': u'A',
                },
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 147,
                    'EXPERIMENT_SAMPLE_ID': None,
                    'EXPERIMENT_SAMPLE_LABEL': None,
                    'EXPERIMENT_SAMPLE_UID': None,
                    'VARIANT_SET_LABEL': ['green'],
                    'VARIANT_SET_UID': ['88349a95'],
                    'ALT': u'C',
                    'REF': u'A',
                }
        ]

        modified_obj_list = _modify_obj_list_for_variant_set_display(
                initial_obj_list)

        # Rows preserved.
        self.assertEqual(2, len(modified_obj_list))


    def test_modify_obj_list_for_variant_set_display__catch_all(self):
        """Test that we preserve the rows for different alts.
        """
        GREEN_VS_LABEL = 'green'
        GREEN_VS_UID = '88349a95'
        initial_obj_list = [
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_UID': None,
                    'VARIANT_SET_LABEL': [GREEN_VS_LABEL],
                    'VARIANT_SET_UID': [GREEN_VS_UID],
                },
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_UID': 'abc',
                    'VARIANT_SET_LABEL': [None],
                    'VARIANT_SET_UID': [None],
                },
                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_UID': 'def',
                    'VARIANT_SET_LABEL': [None],
                    'VARIANT_SET_UID': [None],
                },

                {
                    'ID': 145,
                    'UID': u'7ad36beb',
                    'VA_ID': 146,
                    'EXPERIMENT_SAMPLE_UID': 'ghi',
                    'VARIANT_SET_LABEL': [None],
                    'VARIANT_SET_UID': [None],
                },
        ]

        modified_obj_list = _modify_obj_list_for_variant_set_display(
                initial_obj_list)

        # Rows preserved.
        self.assertEqual(3, len(modified_obj_list))

        for obj in modified_obj_list:
            variant_set_label_list = obj[MELTED_SCHEMA_KEY__VS_LABEL]
            self.assertEqual(1, len(variant_set_label_list))
            self.assertIsNone(variant_set_label_list[0])

            variant_set_uid_list = obj[MELTED_SCHEMA_KEY__VS_UID]
            self.assertEqual(1, len(variant_set_uid_list))
            self.assertIsNone(variant_set_uid_list[0])

            variant_set_label_list = obj[ALL_VS_LABEL_KEY]
            self.assertEqual(1, len(variant_set_label_list))
            self.assertEqual(GREEN_VS_LABEL, variant_set_label_list[0])

            variant_set_uid_list = obj[ALL_VS_UID_KEY]
            self.assertEqual(1, len(variant_set_uid_list))
            self.assertEqual(GREEN_VS_UID, variant_set_uid_list[0])
