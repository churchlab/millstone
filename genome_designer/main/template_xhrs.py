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
from main.models import Project


# Controls ====================================================================
# These are request for model-specific dropdown buttons that go inside the
# dataTable sDom. These templates also includes any matching modals.

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
        'ref_genome': ref_genome,
        'is_melted': True
    }

    context['table_id'] = request.GET.get('tableId', 'sample-list-datatable')

    return HttpResponse(
            render_to_string('controls/variant_filter_controls.html', context))


def reference_genome_list_controls(request):
    """Returns the Reference Genome List control box.
    """
    project_uid = request.GET.get('projectUid')
    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    csrf = request.GET.get('csrf')

    context = {
        'csrf_token': csrf,
        'project': project,
    }

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context['table_id'] = request.GET.get('tableId',
            'reference-genome-list-datatable')

    return HttpResponse(
            render_to_string('controls/reference_genome_list_controls.html',
                    context))


def sample_list_controls(request):
    """Returns the Sample List control box.
    """
    project_uid = request.GET.get('projectUid')
    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    csrf = request.GET.get('csrf')

    context = {
        'csrf_token': csrf,
        'project': project,
    }

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context['table_id'] = request.GET.get('tableId', 'sample-list-datatable')

    return HttpResponse(
            render_to_string('controls/sample_list_controls.html',
            context))


def alignment_list_controls(request):
    """Returns the Alignment List control box.
        Requires no model.
    """
    project_uid = request.GET.get('projectUid')
    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    context = {
        'project': project
    }

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context['table_id'] = request.GET.get('tableId', 'sample-list-datatable')

    return HttpResponse(
            render_to_string('controls/alignment_list_controls.html',
            context))


def variant_set_list_controls(request):
    """Returns the Variant Set List control box.
        Requires no model.
    """
    context = {}

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context['table_id'] = request.GET.get('tableId', 'sample-list-datatable')

    return HttpResponse(
            render_to_string('controls/variant_set_list_controls.html',
            context))
