"""
Methods that handle Ajax requests from the frontend.

This module was created in response to views.py getting quite big, and so a
reasonable separation point is to separate page actions from Ajax actions.
"""

import csv

from django.contrib.auth.decorators import login_required
from django.http import Http404
from django.http import HttpResponse
from django.shortcuts import get_object_or_404

from main.models import VariantSet
from scripts.variant_filter import get_variants_that_pass_filter


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

