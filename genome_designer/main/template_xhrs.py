"""
Methods that return template html.

We use these for lack of integration of a templating solution that is supported
in both python and javascript.  Front-end views request the template and get
the html string as the response.
"""

from django.conf import settings
from django.http import HttpResponse
from django.shortcuts import get_list_or_404
from django.shortcuts import get_object_or_404
from django.template import RequestContext
from django.template.loader import render_to_string

from main.models import AlignmentGroup
from main.models import ReferenceGenome
from main.models import Project
from main.models import SavedVariantFilterQuery


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

    saved_query_list = SavedVariantFilterQuery.objects.filter(
            owner=request.user.get_profile())

    example_saved_query_list = [
        'GT_TYPE = 2',
        'INFO_EFF_IMPACT = HIGH'
    ]

    context = {
        'ref_genome': ref_genome,
        'is_melted': True,
        'saved_query_list': saved_query_list,
        'example_saved_query_list': example_saved_query_list
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
        'project': project,
    }

    # If there are no de novo assemblies, pass flag to context to hide the
    # toggle de novo assemblies button
    project_has_de_novo_assemblies = False
    for rg in ReferenceGenome.objects.filter(project=project):
        if rg.metadata.get('is_from_de_novo_assembly', False):
            project_has_de_novo_assemblies = True
            break
    context['project_has_de_novo_assemblies'] = project_has_de_novo_assemblies

    context['show_de_novo'] = int(request.GET.get('showDeNovo'))

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context['table_id'] = request.GET.get('tableId',
            'reference-genome-list-datatable')

    return HttpResponse(
            render_to_string('controls/reference_genome_list_controls.html',
                    context))


def contig_list_controls(request):
    """Returns the Contig List control box.
    """

    # TODO: Find out if below line is important
    # csrf = request.GET.get('csrf')

    alignment_group = get_object_or_404(
            AlignmentGroup,
            uid=request.GET.get('alignmentGroupUid'))

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context = {
            'table_id': request.GET.get('tableId',
                    'reference-genome-list-datatable'),
            'alignment_group_uid': alignment_group.uid,
            'samples_uid_tuples': [
                    (esta.experiment_sample.label, esta.uid) for esta in
                    alignment_group.experimentsampletoalignment_set.all()]
    }

    return HttpResponse(
            render_to_string('controls/contig_list_controls.html', context))


def sample_list_controls(request):
    """Returns the Sample List control box.
    """
    project_uid = request.GET.get('projectUid')
    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    csrf = request.GET.get('csrf')

    context = {
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

def alignment_controls(request):
    """Returns the single Alignment control box.
        Requires no model.
    """
    project_uid = request.GET.get('projectUid')
    alignment_group_uid = request.GET.get('alignmentGroupUid')

    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    alignment_group = get_object_or_404(AlignmentGroup,
            reference_genome__project__owner=request.user.get_profile(),
            uid=alignment_group_uid)

    context = {
        'FLAG__GENOME_FINISHING_ENABLED': settings.FLAG__GENOME_FINISHING_ENABLED,
        'project': project,
        'table_id': request.GET.get('tableId'),
        'alignment_group': alignment_group
    }

    return HttpResponse(
            render_to_string('controls/alignment_controls.html',
            context))


def variant_set_list_controls(request):
    """Returns the Variant Set List control box.
    """
    # Validate the request by making sure this is a valid project for the user.
    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=request.GET.get('projectUid'))

    context = RequestContext(request)

    # Get the list of ReferenceGenomes for the project for the VariantSet
    # creation modal.
    context['ref_genome_list'] = ReferenceGenome.objects.filter(project=project)

    # If the request passed a tableId, then give it to Django to decorate the
    # controls.
    context['table_id'] = request.GET.get('tableId', 'sample-list-datatable')

    return HttpResponse(
            render_to_string('controls/variant_set_list_controls.html',
            context))


def create_new_empty_variant_set(request):
    """Creating a new empty variant set
    """
    ref_genome_uid = request.GET.get('refGenomeUid')

    ref_genome_list = [get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)]

    context = {
        'ref_genome_list': ref_genome_list
    }

    return HttpResponse(render_to_string(
        'controls/create_empty_variant_set_modal.html', context))

