"""
Tests for model_views.py.
"""

from django.test import TestCase

from main.model_views import _modify_obj_list_for_variant_set_display


class TestModelViews(TestCase):

    def test_modify_obj_list_for_variant_set_display__different_alts(self):
        """Test that we preserve the rows for different alts.
        """
        initial_obj_list = [
                {'VA_ID': 146, 'VA_DATA': {}, 'VE_ID': None, 'UID': u'7ad36beb', 'EXPERIMENT_SAMPLE_LABEL': None, 'EXPERIMENT_SAMPLE_UID': None, 'VARIANT_SET_LABEL': ['green'], 'VARIANT_SET_UID': ['88349a95'], 'VCCD_ID': None, 'EXPERIMENT_SAMPLE_ID': None, 'POSITION': 1313L, 'ALT': u'T', 'REF': u'A', 'ID': 145, 'CHROMOSOME': u'c'},
                {'VA_ID': 147, 'VA_DATA': {}, 'VE_ID': None, 'UID': u'7ad36beb', 'EXPERIMENT_SAMPLE_LABEL': None, 'EXPERIMENT_SAMPLE_UID': None, 'VARIANT_SET_LABEL': ['green'], 'VARIANT_SET_UID': ['88349a95'], 'VCCD_ID': None, 'EXPERIMENT_SAMPLE_ID': None, 'POSITION': 1313L, 'ALT': u'C', 'REF': u'A', 'ID': 145, 'CHROMOSOME': u'c'}
        ]
        modified_obj_list = _modify_obj_list_for_variant_set_display(
                initial_obj_list)
        self.assertEqual(2, len(modified_obj_list))
