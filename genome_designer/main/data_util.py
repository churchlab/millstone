"""
Common methods for getting data from the backend.

These methods are intended to be used by both views.py, which should define
only pages, and xhr_handlers.py, which are intended to respond to AJAX
requests.
"""

from main.model_views import CastVariantView
from main.model_views import MeltedVariantView
from variants.variant_filter import get_variants_that_pass_filter


def lookup_variants(reference_genome, combined_filter_string, is_melted):
    """Lookup the Variants that match the filter specified in the params.

    Returns:
        List of CastVariantView or MeltedVariantView objects.
    """
    # Apply the filters.
    filter_result = get_variants_that_pass_filter(
            combined_filter_string, reference_genome)
    variant_list = filter_result.variant_set
    variant_id_to_metadata_dict = filter_result.variant_id_to_metadata_dict

    # Convert to appropriate view objects.
    if is_melted:
        melted_variant_list = []
        for variant in variant_list:
            melted_variant_list.extend(
                    MeltedVariantView.variant_as_melted_list(variant,
                            variant_id_to_metadata_dict))
        return melted_variant_list
    else:
        return [CastVariantView.variant_as_cast_view(variant,
                variant_id_to_metadata_dict) for variant in variant_list]
