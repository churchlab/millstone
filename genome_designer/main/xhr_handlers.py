"""
Methods that handle Ajax requests from the frontend.

This module was created in response to views.py getting quite big, and so a
reasonable separation point is to separate page actions from Ajax actions.
"""

import csv
import json

from django.contrib.auth.decorators import login_required
from django.http import Http404
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_http_methods

from main.adapters import adapt_model_to_frontend
from main.data_util import lookup_variants
from main.data_util import VARIANT_FILTER_STRING_KEY
from main.models import Project
from main.models import ReferenceGenome
from main.models import VariantSet
from scripts.variant_filter import get_variants_that_pass_filter
from scripts.variant_sets import add_or_remove_variants_from_set


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

    # Get inputs to perform the query for Variants data.
    if VARIANT_FILTER_STRING_KEY in request.GET:
        manual_filter_string = request.GET.get(VARIANT_FILTER_STRING_KEY)
    else:
        manual_filter_string = ''
    # TODO: Combine with saved filter string.
    combined_filter_string = manual_filter_string
    is_melted = request.GET.get('melt', 0) == '1'

    response_data = {
        'variant_list_json': lookup_variants(reference_genome,
                combined_filter_string, is_melted),
        'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                {'reference_genome__project': project}),
    }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


@login_required
@require_http_methods(['POST'])
def modify_variant_in_set_membership(request):
    """Action that handles modifying the membership of a Variant in a
    VariantSet.
    """
    request_data = json.loads(request.body)

    # Parse the data from the request body.

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
    # project = get_object_or_404(Project, owner=request.user.get_profile(),
    #         uid=project_uid)
    project = Project.objects.get(owner=request.user.get_profile(),
            uid=project_uid)
    reference_genome = ReferenceGenome.objects.get(project=project,
            uid=ref_genome_uid)


    context = {}

    # Add or remove the variants to the set, as variantSetAction.
    context.update(add_or_remove_variants_from_set(
        request_data.get('variantUidList'),
        request_data.get('variantSetAction'),
        request_data.get('variantSetUid'),
    ))

    # TODO: Pass these from the frontend.
    combined_filter_string = ''
    is_melted = 0

    context['variant_list_json'] = lookup_variants(reference_genome,
        combined_filter_string, is_melted)
    context['variant_set_list_json'] = adapt_model_to_frontend(VariantSet,
            {'reference_genome__project': project})

    return HttpResponse(
            json.dumps(context),
            content_type='application/json')
