"""
Common methods for getting data from the backend.

These methods are intended to be used by both views.py, which should define
only pages, and xhr_handlers.py, which are intended to respond to AJAX
requests.

This module interacts closely with the ModelViews in model_views.py.
"""

from collections import defaultdict

from django.db import connection

from main.model_views import CastVariantView
from main.model_views import MeltedVariantView
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantEvidence
from variants.common import dictfetchall
from variants.materialized_variant_filter import get_variants_that_pass_filter


class LookupVariantsResult(object):
    """Result of a call to lookup_variants.

    Attributes:
        result_list: List of cast or melted Variant objects.
        num_total_variants: Total number of variants that match query.
            For pagination.
    """
    def __init__(self, result_list, num_total_variants):
        self.result_list = result_list
        self.num_total_variants = num_total_variants


def lookup_variants(reference_genome, combined_filter_string, is_melted,
        pagination_start, pagination_len):
    """Manages the end-to-end flow of looking up Variants that match the
    given filter.

    This function delegates to the variant_filter module to get the list of
    variants matching the filter. Then, this function takes those results
    and handles casting them to appropriate view-type objects (e.g. Melted vs
    Cast).

    Returns:
        LookupVariantsResult object.
    """
    # First get the Variants that pass the filter.
    filter_eval_result = get_variants_that_pass_filter(combined_filter_string,
            reference_genome)
    result_list = list(filter_eval_result.variant_set)

    # If this is a melted view, return results as they are.
    if is_melted:
        # TODO: Handle pagination.
        page_results = result_list[pagination_start :
                pagination_start + pagination_len]
        num_total_variants = 1000000
        return LookupVariantsResult(page_results, num_total_variants)

    # Otherwise, we need to Cast the results.
    page_results = cast_joined_variant_objects(result_list)
    page_results = page_results[pagination_start :
            pagination_start + pagination_len]
    num_total_variants = 1000000
    return LookupVariantsResult(page_results, num_total_variants)


def cast_joined_variant_objects(melted_variant_list):
    """Converts the list of melted variants into a cast representation.

    This means returning one row per variant, compressing other columns
    into an aggregate representation. For example, the 'experiment_sample_uid'
    column becomes the 'total_samples'.
    """
    cast_obj_list = []

    # First, we build a structure from variant id to list of result rows.
    variant_id_to_result_row = defaultdict(list)
    for result in melted_variant_list:
        variant_id_to_result_row[result['id']].append(result)

    for variant_id, result_row_list in variant_id_to_result_row.iteritems():
        assert len(result_row_list), "Not expected. Debug."
        position = result_row_list[0]['position']
        ref = result_row_list[0]['ref']
        total_samples = len(result_row_list)
        cast_obj_list.append({
            'id': variant_id,
            'position': position,
            'ref': ref,
            'alt': 'TODO',
            'total_samples': total_samples
        })

    return cast_obj_list
