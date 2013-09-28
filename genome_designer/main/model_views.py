"""
Classes that describe how a particular model should be viewed.
"""

from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from scripts.dynamic_snp_filter_key_map import MAP_KEY__COMMON_DATA
from scripts.dynamic_snp_filter_key_map import MAP_KEY__EVIDENCE


class BaseModelView(object):
    """Common methods for model views.

    NOTE: May want to look into more of a mixin model instead.
    """
    pass


class CastVariantView(BaseModelView):
    """View of Variant that does not expose the individual
    VariantToExperimentSample relationships.
    """
    pass


class MeltedVariantView(BaseModelView):
    """View of a Variant when showing one row per VariantToExperimentSample.
    """

    def __init__(self, variant, variant_caller_common_data=None, variant_evidence=None):
        self.variant = variant
        self.variant_caller_common_data = variant_caller_common_data
        self.variant_evidence = variant_evidence

    def custom_getattr(self, attr):
        delegate_order = [
                self.variant,
                self.variant_caller_common_data,
                self.variant_evidence
        ]
        for delegate in delegate_order:
            if delegate is not None and hasattr(delegate, attr):
                return getattr(delegate, attr)
        # Default.
        return 'undefined'

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.

        Args:
            field_mask: A dictionary from field name to True/False that defines
                whether to include that key. Otherwise, the default is
                returned.
        """
        if 'visible_key_names' in kwargs:
            variant_key_map = kwargs['variant_key_map']
            visible_key_names = kwargs['visible_key_names']

            additional_common_data_field_list = [field for field in
                    visible_key_names if field in
                    variant_key_map[MAP_KEY__COMMON_DATA]]

            additional_evidence_field_list = [field for field in visible_key_names
                    if field in variant_key_map[MAP_KEY__EVIDENCE]]
        else:
            additional_common_data_field_list = []
            additional_evidence_field_list = []

        # Sort the visible keys into appropriate buckets based on the map.
        return (Variant.get_field_order() +
                VariantCallerCommonData.get_field_order(
                        additional_field_list=additional_common_data_field_list) +
                VariantEvidence.get_field_order(
                        additional_field_list=additional_evidence_field_list))
