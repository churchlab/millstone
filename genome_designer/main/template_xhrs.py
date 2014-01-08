"""
Methods that return template html.

We use these for lack of integration of a templating solution that is supported
in both python and javascript.  Front-end views request the template and get
the html string as the response.
"""

from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.template.loader import render_to_string

from main.models import ReferenceGenome


def variant_filter_controls(request):
    """Returns the Variant filter control box.
    """
    ref_genome_uid = request.GET.get('refGenomeUid')
    ref_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    csrf = request.GET.get('csrf')

    context = {
        'csrf': csrf,
        'ref_genome': ref_genome
    }

    return HttpResponse(
            render_to_string('snp_filter_control.html', context))


def variant_set_controls(request):
    """Returns the VariantSets control box.
    """
    return HttpResponse(render_to_string('variant_set_controls.html'))
