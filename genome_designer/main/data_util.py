"""
Common methods for getting data from the backend.

These methods are intended to be used by both views.py, which should define
only pages, and xhr_handlers.py, which are intended to respond to AJAX
requests.

This module interacts closely with the ModelViews in model_views.py.
"""

from collections import defaultdict
from itertools import groupby

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


def lookup_variants(query_args, reference_genome):
    """Manages the end-to-end flow of looking up Variants that match the
    given filter.

    This function delegates to the variant_filter module to get the list of
    variants matching the filter, and in cast or melted form depending on the
    is_melted flag.

    Returns:
        LookupVariantsResult object.
    """
    # Get the Variants that pass the filter.
    page_results = get_variants_that_pass_filter(query_args, reference_genome)

    # Now get all results that pass the filter (remove limit clause)
    query_args['count_only'] = True
    query_args['pagination_start'] = 0
    query_args['pagination_len'] = -1  # for no limit
    num_total_variants = get_variants_that_pass_filter(query_args, reference_genome)[0]['count']

    # Maybe cast the results.
    if not query_args['is_melted']:
        page_results = format_cast_objects(page_results)

    return LookupVariantsResult(page_results, num_total_variants)


def format_cast_objects(page_results):
    for page_result in page_results:
        # If there is an empty row (no ExperimentSample associated),
        # then all the counts will be off by one, so we need to decrement
        # them.
        if None in page_result['experiment_sample_uid']:
            maybe_dec = 1
        else:
            maybe_dec = 0
        assert maybe_dec >= 0, "maybe_dec should be positive"

        alts = page_result['alt']

        # Append ref with count of # variants without alt
        page_result['ref'] += ' (%d)' % alts.count(None)

        # List frequency of each alt
        processed_alts = sorted(filter(lambda alt: alt, alts))
        page_result['alt'] = ' | '.join(['%s (%d)' %
            (val, len(list(group)) - maybe_dec)
            for val, group in groupby(processed_alts)])

        # Combine variant sets, with frequencies.
        variant_set_label_list = page_result['variant_set_label']
        processed_variant_set_label_list = sorted(filter(
                lambda vs_label: vs_label, variant_set_label_list))
        page_result['variant_sets'] = ' | '.join(['%s (%d)' %
            (val, len(list(group)) - maybe_dec) for val, group in groupby(
                    processed_variant_set_label_list)])

        # Add additional information specific to cast view
        page_result['total_samples'] = len(alts) - maybe_dec
    return page_results


def cast_object_list_field_as_bucket_string(cast_object_dict_list, field):
    """Converts a Cast object's field as a list into a string where values
    have been bucketed by unique type, with counts in parens.

    TODO: This doesn't really make sense for fields that take on continuous
    number values. Figure out what to do with these kinds of fields.
    """
    # First, extract the relevant fields.
    value_list = []
    for cast_object_dict in cast_object_dict_list:
        if cast_object_dict is None:
            value_list.append('NONE')
        else:
            value_list.append(cast_object_dict.get(field, ''))

    # Count.
    buckets = defaultdict(lambda: 0)
    for val in value_list:
        buckets[str(val)] += 1

    # Create string.
    return (' | '.join(['%s (%d)' % (key, count)
            for key, count in buckets.iteritems()]))
