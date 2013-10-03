"""
Classes that describe how a particular model should be viewed.
"""

from django.db.models.query import QuerySet

from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from scripts.dynamic_snp_filter_key_map import MAP_KEY__COMMON_DATA
from scripts.dynamic_snp_filter_key_map import MAP_KEY__EVIDENCE


class BaseVariantView(object):
    """Common methods for model views.
    """

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


class CastVariantView(BaseVariantView):
    """View of Variant that does not expose the individual
    VariantToExperimentSample relationships.
    """
    def __init__(self, variant):
        self.variant = variant
        self.variant_caller_common_data_list = (
                variant.variantcallercommondata_set.all())
        self.variant_evidence_list = VariantEvidence.objects.filter(
                variant_caller_common_data__variant=variant)

    def custom_getattr(self, attr):
        """For an attribute of a Variant, returns what you would expect.

        For many-to-one relations, e.g. VariantCallerCommonData and
        VariantEvidence, returns a '|'-separated list of values for
        that attribute.

        Returns:
            A string representing the view for the attribute.
        """
        delegate_order = [
                self.variant,
                self.variant_caller_common_data_list,
                self.variant_evidence_list
        ]
        for delegate in delegate_order:
            result_list = []

            # Convert to list for next step.
            if isinstance(delegate, QuerySet):
                delegate_list = delegate
            else:
                delegate_list = [delegate]

            # Iterate through the delegates.
            for delegate in delegate_list:

                if hasattr(delegate, attr):
                    result_list.append(getattr(delegate, attr))
                else:
                    try: result_list.append(delegate.as_dict()[attr])
                    except: pass

            if len(result_list) == 1:
                # If one object, return it.
                return result_list[0]
            elif len(result_list) > 1:
                # Else if more than one, return a '|'-separated string.
                return ' | '.join([str(res) for res in result_list])

        # Default.
        return 'undefined'

    @classmethod
    def variant_as_cast_view(clazz, variant_obj,
            variant_id_to_metadata_dict=None):
        """Factory method returns a cast view for the given variant.

        Args:
            variant_obj: The Variant object to melt.
            variant_id_to_metadata_dict: See TODO.

        Returns:
            A CastVariantView instance.
        """
        # TODO: Do we need to modify the view based on passing samples in the
        # metadata.
        return CastVariantView(variant_obj)


class MeltedVariantView(BaseVariantView):
    """View of a Variant when showing one row per VariantToExperimentSample.
    """

    def __init__(self, variant, variant_caller_common_data=None,
            variant_evidence=None):
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
            if delegate is not None:
                if hasattr(delegate, attr):
                    return getattr(delegate, attr)
                else:
                    try:
                        return delegate.as_dict()[attr]
                    except:
                        pass
        # Default.
        return '---'

    @classmethod
    def variant_as_melted_list(clazz, variant_obj,
            variant_id_to_metadata_dict=None):
        """Factory method that melts the variant into a list of objects,
        one per common data per sample.

        Args:
            variant_obj: The Variant object to melt.
            variant_id_to_metadata_dict: Restrict melted entities to those
                related to these samples.

        Returns:
            List of MeltedVariantView objects corresponding to the Variant.
            This list will always be at least length 1.
        """
        melted_list = []

        all_common_data = variant_obj.variantcallercommondata_set.all()

        # We guarantee returning at least one row.
        if len(all_common_data) == 0:
            melted_list.append(MeltedVariantView(variant_obj, None, None))
            return melted_list

        # Otherwise iterate through
        for common_data_obj in all_common_data:
            # Return one row not associated with any sample.
            melted_list.append(MeltedVariantView(variant_obj, common_data_obj, None))

            variant_evidence_list = common_data_obj.variantevidence_set.all()
            if len(variant_evidence_list) == 0:
                continue

            if variant_id_to_metadata_dict is not None:
                passing_sample_ids = variant_id_to_metadata_dict[variant_obj.id].get(
                        'passing_sample_ids', set())
            else:
                passing_sample_ids = None

            for variant_evidence_obj in variant_evidence_list:
                if passing_sample_ids is not None:
                    sample_id = variant_evidence_obj.experiment_sample.id
                    if not sample_id in passing_sample_ids:
                        continue
                melted_list.append(MeltedVariantView(variant_obj, common_data_obj,
                        variant_evidence_obj))

        return melted_list
