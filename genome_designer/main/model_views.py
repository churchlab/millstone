"""
Classes that describe how a particular model should be viewed.
"""

from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantEvidence


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

    def __init__(self, variant, variant_caller_common_data, variant_evidence):
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
            if hasattr(delegate, attr):
                return getattr(delegate, attr)
        raise AttributeError("MeltedVariantView object has not attribute %s" % attr)

    @classmethod
    def get_field_order(clazz):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return (Variant.get_field_order() +
                VariantCallerCommonData.get_field_order() +
                VariantEvidence.get_field_order())
