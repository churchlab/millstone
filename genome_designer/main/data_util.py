"""
Common methods for getting data from the backend.

These methods are intended to be used by both views.py, which should define
only pages, and xhr_handlers.py, which are intended to respond to AJAX
requests.
"""

from main.adapters import adapt_model_or_modelview_list_to_frontend
from main.melt_util import variant_as_melted_list
from scripts.variant_filter import get_variants_that_pass_filter


# Key in the GET params containing the string for filtering the variants.
VARIANT_FILTER_STRING_KEY = 'variantFilterString'


def lookup_variants(reference_genome, combined_filter_string, is_melted):
    """Lookup the Variants that match the filter specified in the params.
    """
    # Apply the filters.
    filter_result = get_variants_that_pass_filter(
            combined_filter_string, reference_genome)
    variant_list = filter_result.variant_set
    variant_id_to_metadata_dict = filter_result.variant_id_to_metadata_dict

    # Determine whether the melted or cast version is requested.
    if is_melted:
        melted_variant_list = []
        for variant in variant_list:
            melted_variant_list.extend(variant_as_melted_list(variant,
                    variant_id_to_metadata_dict))
        variant_list = melted_variant_list

    # TODO: I'm confused where we have sets and where we have lists.
    variant_list = list(variant_list)

    return adapt_model_or_modelview_list_to_frontend(variant_list)
