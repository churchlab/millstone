"""
Methods that handle Ajax requests from the frontend.

This module was created in response to views.py getting quite big, and so a
reasonable separation point is to separate page actions from Ajax actions.
"""

import copy
import csv
import json

from django.contrib.auth.decorators import login_required
from django.http import Http404
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_http_methods

from main.adapters import adapt_model_or_modelview_list_to_frontend
from main.adapters import adapt_model_to_frontend
from main.data_util import lookup_variants
from main.models import Project
from main.models import ReferenceGenome
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from scripts.dynamic_snp_filter_key_map import MAP_KEY__COMMON_DATA
from scripts.dynamic_snp_filter_key_map import MAP_KEY__EVIDENCE
from variants.common import extract_filter_keys
from variants.variant_filter import get_variants_that_pass_filter
from variants.variant_sets import add_or_remove_variants_from_set


@login_required
def export_variant_set_as_csv(request):
    """Returns a response to download a file for all of
    the samples contained in the variant set.
    """
    variant_set_uid = request.POST.get('variant_set_uid', None)
    if not variant_set_uid:
        raise Http404

    # Make sure that the variant set is in a Project belonging to the user.
    variant_set = get_object_or_404(VariantSet,
            uid=variant_set_uid,
            reference_genome__project__owner=request.user.get_profile())

    # Get the variants passing the filter.
    filter_string = 'IN_SET(%s)' % variant_set_uid
    filter_result = get_variants_that_pass_filter(
            filter_string, variant_set.reference_genome)
    variant_list = filter_result.variant_set

    # Export data as csv.
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="variant_set.csv"'
    csv_field_names = [
        'position',
        'ref',
        'alt'
    ]
    writer = csv.DictWriter(response, csv_field_names)
    writer.writeheader()
    for variant in variant_list:
        writer.writerow({
            'position': variant.position,
            'ref': variant.ref_value,
            'alt': variant.alt_value,
        })
    return response


# Key in the GET params containing the string for filtering the variants.
VARIANT_FILTER_STRING_KEY = 'variantFilterString'

@login_required
def get_variant_list(request):
    """Returns a list of Variants, filtered by any filter parameters contained
    in the request.
    """
    # Parse the GET params.
    ref_genome_uid = request.GET.get('refGenomeUid')
    project_uid = request.GET.get('projectUid')

    # Get the project and verify that the requesting user has the
    # right permissions.
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    reference_genome = ReferenceGenome.objects.get(project=project,
            uid=ref_genome_uid)

    # Pagination.
    pagination_start = int(request.GET.get('iDisplayStart', 0))
    pagination_len = int(request.GET.get('iDisplayLength', 100))

    # Get inputs to perform the query for Variants data.
    if VARIANT_FILTER_STRING_KEY in request.GET:
        manual_filter_string = request.GET.get(VARIANT_FILTER_STRING_KEY)
    else:
        manual_filter_string = ''
    # TODO: Combine with saved filter string.
    combined_filter_string = manual_filter_string
    is_melted = request.GET.get('melt', 0) == '1'

    # Determine the visible keys.
    visible_key_names = _determine_visible_field_names(request,
            combined_filter_string, reference_genome)

    # Get the list of Variants (or melted representation) to display.
    lookup_variant_result = lookup_variants(reference_genome, combined_filter_string,
            is_melted, pagination_start, pagination_len)
    variant_list = lookup_variant_result.result_list
    num_total_variants = lookup_variant_result.num_total_variants
    variant_list_json = adapt_model_or_modelview_list_to_frontend(variant_list,
            variant_key_map=reference_genome.variant_key_map,
            visible_key_names=visible_key_names)

    # Grab the VariantSet data.
    variant_set_list = VariantSet.objects.filter(
            reference_genome=reference_genome)

    # Query the keys valid for ReferenceGenome, and mark the ones that
    # will be displayed so that the checkmarks are pre-filled in case
    # the user wishes to change these.
    variant_key_map_with_active_fields_marked = copy.deepcopy(
            reference_genome.variant_key_map)
    _mark_active_keys_in_variant_key_map(
            variant_key_map_with_active_fields_marked)

    response_data = {
        'variant_list_json': variant_list_json,
        'num_total_variants': num_total_variants,
        'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                obj_list=variant_set_list),
        'variant_key_filter_map_json': json.dumps(
                variant_key_map_with_active_fields_marked)
    }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


def _determine_visible_field_names(request, filter_string, ref_genome):
    """Determine which fields to show.
    """
    # Get visible keys explicitly marked in the UI by the user.
    if 'visibleKeyNames' in request.GET:
        visible_key_names = json.loads(request.GET.get('visibleKeyNames'))
    else:
        visible_key_names = []

    # Also show keys in the filter string.
    fields_from_filter_string = extract_filter_keys(filter_string, ref_genome)

    return list(set(visible_key_names) | set(fields_from_filter_string))


def _mark_active_keys_in_variant_key_map(variant_key_map):
    """Mutates variant_key_map to mark fields that should be active based
    on model class defaults.
    """
    # In the current implementation, we mark the fields that are included
    # in the relevant models' get_field_order() method.

    def _update_model_class_key_map(model_class, variant_key_submap):
        """Helper method."""
        default_keys = [el['field'] for el in model_class.get_field_order()]
        for key in default_keys:
            if key in variant_key_submap:
                variant_key_submap[key]['checked'] = True

    _update_model_class_key_map(VariantCallerCommonData,
            variant_key_map[MAP_KEY__COMMON_DATA])

    _update_model_class_key_map(VariantEvidence,
            variant_key_map[MAP_KEY__EVIDENCE])


@login_required
@require_http_methods(['POST'])
def modify_variant_in_set_membership(request):
    """Action that handles modifying the membership of a Variant in a
    VariantSet.
    """
    request_data = json.loads(request.body)

    # Make sure the required keys are present.
    REQUIRED_KEYS = [
            'variantUidList',
            'variantSetAction',
            'variantSetUid']

    # Validate the request.
    if not all(key in request_data for key in REQUIRED_KEYS):
        return HttpResponseBadRequest("Invalid request. Missing keys.")

    ref_genome_uid = request_data.get('refGenomeUid')
    project_uid = request_data.get('projectUid')

    # Get the project and verify that the requesting user has the
    # right permissions.
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    reference_genome = get_object_or_404(ReferenceGenome, project=project,
            uid=ref_genome_uid)

    context = {}

    # Add or remove the variants to the set, as variantSetAction.
    context.update(add_or_remove_variants_from_set(
        request_data.get('variantUidList'),
        request_data.get('variantSetAction'),
        request_data.get('variantSetUid'),
    ))

    # TODO: Proper response and error handling.

    return HttpResponse('ok')
