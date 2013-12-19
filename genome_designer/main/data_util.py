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
from variants.variant_filter import get_variants_that_pass_filter


class RequestScopedVariantDataCache(object):
    """Data cached scoped to the request.

    This object is used to make bulk data calls for all data that will be
    needed to respond to a request for Variants.
    """

    def __init__(self):
        self.v_id_to_alt_map = {}
        self.v_id_to_vccd_map = {}
        self.v_id_to_ve_map= {}


    def populate(self, variant_list, project):
        """Populates the cache.
        """
        self.v_id_to_alt_map = defaultdict(list)
        self.v_id_to_vccd_map = defaultdict(list)
        self.v_id_to_ve_map = defaultdict(list)

        # Perform 2 different queries for simplicity. I welcome any more
        # optimal solutions.
        # TODO: Work on this.

        variant_id_list = [variant.id for variant in variant_list]

        # Get all Variant, VariantCallerCommon,  and VariantAlternate data.
        all_variant_data_list = Variant.objects.filter(id__in=variant_id_list).\
                prefetch_related('variantalternate_set').\
                prefetch_related('variantcallercommondata_set')

        # Get VariantEvidence and VariantCallerCommonData.
        all_ve_data_list = VariantEvidence.objects.filter(
                variant_caller_common_data__variant_id__in=variant_list).\
            select_related('variant_caller_common_data')

        # HACK: Set the sample_uid for each ve manually. Is there a better
        # way to do this? Using select_related() above is too slow.
        # See next HACK note below.
        sample_id_to_uid = {}
        for sample in ExperimentSample.objects.filter(project=project):
            sample_id_to_uid[sample.id] = sample.uid

        # Sort the data into the cache object's data structures.
        for v in all_variant_data_list:
            self.v_id_to_alt_map[v.id] = v.variantalternate_set.all()
            self.v_id_to_vccd_map[v.id] = v.variantcallercommondata_set.all()

        for ve in all_ve_data_list:
            variant_id = ve.variant_caller_common_data.variant_id
            self.v_id_to_ve_map[variant_id].append(ve)

            # HACK: Continued from above.
            ve.manually_cached_data['sample_uid'] = (
                    sample_id_to_uid[ve.experiment_sample_id])


    def get_variant_alternate_list(self, variant):
        return self.v_id_to_alt_map[variant.id]


    def get_variant_caller_common_data_list(self, variant):
        return self.v_id_to_vccd_map[variant.id]


    def get_variant_evidence_list(self, variant):
        return self.v_id_to_ve_map[variant.id]


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
    """Lookup the Variants that match the filter specified in the params.

    Returns:
        LookupVariantsResult object.
    """
    # First get the Variants that pass the filter.

    cursor = connection.cursor()
    # We can do a rough limit on how many results.
    # Estimate the total number of results as the number of variants
    # cross the number of samples.
    # NOTE: This should be an over-calculation, I believe.
    # TODO: Figure this out precisely.
    num_samples = ExperimentSample.objects.filter(
        project=reference_genome.project).count()
    num_pages = pagination_start / pagination_len
    approx_num_max_results = (num_pages + 1) * pagination_len * num_samples

    sql_statement = (
            'SELECT id, position, experiment_sample_uid '
            'FROM materialized_melted_variant '
            'LIMIT %d'
            % (approx_num_max_results,)
    )
    cursor.execute(sql_statement)
    result_list = dictfetchall(cursor)

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
        total_samples = len(result_row_list)
        cast_obj_list.append({
            'id': variant_id,
            'position': position,
            'total_samples': total_samples
        })

    return cast_obj_list
