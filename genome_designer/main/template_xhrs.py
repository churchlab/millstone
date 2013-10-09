"""
Methods that return template html.

We use these for lack of integration of a templating solution that is supported
in both python and javascript.  Front-end views request the template and get
the html string as the response.
"""

from django.http import HttpResponse
from django.template.loader import render_to_string

def variant_filter_controls(request):
    """Returns the Variant filter control box.
    """
    return HttpResponse(render_to_string('snp_filter_control.html'))


def variant_set_controls(request):
    """Returns the VariantSets control box.
    """
    return HttpResponse(render_to_string('variant_set_controls.html'))
