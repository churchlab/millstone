"""
Common methods for getting data from the backend.

These methods are intended to be used by both views.py, which should define
only pages, and xhr_handlers.py, which are intended to respond to AJAX
requests.
"""

from main.model_views import CastVariantView
from main.model_views import MeltedVariantView
from variants.variant_filter import get_variants_that_pass_filter


class LookupVariantsResult(object):
    """Result of a call to lookup_variants.
    """
    def __init__(self, result_list, num_total_variants):
        self.result_list = result_list
        self.num_total_variants = num_total_variants


def lookup_variants(reference_genome, combined_filter_string, is_melted,
        pagination_start, pagination_len):
    """Lookup the Variants that match the filter specified in the params.

    Returns:
        List of CastVariantView or MeltedVariantView objects.
    """
    # Apply the filters.
    filter_result = get_variants_that_pass_filter(
            combined_filter_string, reference_genome)
    variant_list = list(filter_result.variant_set)
    variant_id_to_metadata_dict = filter_result.variant_id_to_metadata_dict

    # Convert to appropriate view objects.
    if is_melted:
        result_list = []
        for variant in variant_list:
            result_list.extend(
                    MeltedVariantView.variant_as_melted_list(variant,
                            variant_id_to_metadata_dict))

        # Count the results and return just the page we are looking at now.
        num_total_variants = len(result_list)
        page_results = result_list[
            pagination_start:pagination_start + pagination_len]

    else:
        num_total_variants = len(variant_list)

        # Apply pagination limits here before casting since the number is 1:1.
        # TODO: Figure out how to similarly avoid casting all results in the
        # is_melted case above.
        variant_list = variant_list[
                pagination_start:pagination_start + pagination_len]

        # TODO: Avoid casting all results when only subset is returned.
        # Question is how to get the total count without doing the cast.
        page_results = [CastVariantView.variant_as_cast_view(variant,
                variant_id_to_metadata_dict) for variant in variant_list]

    return LookupVariantsResult(page_results, num_total_variants)
