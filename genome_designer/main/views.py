"""
Views of pages.

NOTE: Put new Ajax-only actions into main/xhr_handlers.py.
"""

import json
import os

from django.contrib.auth.decorators import login_required
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.urlresolvers import reverse
from django.http import Http404
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.shortcuts import render

from main.adapters import adapt_model_instance_to_frontend
from main.adapters import adapt_model_or_modelview_list_to_frontend
from main.adapters import adapt_model_to_frontend
from main.forms import ProjectForm
from main.model_views import MeltedVariantView
from main.models import AlignmentGroup
from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.alignment_pipeline import create_alignment_groups_and_start_alignments
from scripts.import_util import import_reference_genome_from_local_file
from scripts.import_util import import_reference_genome_from_ncbi
from scripts.import_util import import_samples_from_targets_file
from scripts.import_util import import_variant_set_from_vcf
from scripts.snp_callers import run_snp_calling_pipeline
import settings

# Tags used to indicate which tab we are on.
TAB_ROOT__DATA = 'DATA'
TAB_ROOT__ANALYZE = 'ANALYZE'


class InputError(Exception):
    """Exception raised for errors in the input.

    I assume this is what DBG wanted in d922023e8a1d4e7757bfa3e68aba3ace89f2a749.

    Attributes:
        msg  -- explanation of the error
    """
    def __init__(self, msg):
        self.msg = msg


def home_view(request):
    """The main landing page.
    """
    context = {}
    return render(request, 'home.html', context)


@login_required
def project_list_view(request):
    """The list of projects.
    """
    # NOTE: The 'project_list' template variable is provided via
    # the custom context processor main.context_processors.common_data.
    context = {}
    return render(request, 'project_list.html', context)


@login_required
def project_create_view(request):
    """View where a user creates a new Project.
    """
    if request.POST:
        project = Project(owner=request.user.get_profile())
        form = ProjectForm(request.POST, instance=project)
        if form.is_valid():
            form.title = form.cleaned_data['title']
            form.save()
            return HttpResponseRedirect(
                    reverse('genome_designer.main.views.project_view',
                            args=(project.uid,)))
    else:
        form = ProjectForm()

    # form is either the form with errors from above or an empty instance.
    context = {
        'form': form,
    }
    return render(request, 'project_create.html', context)


@login_required
def project_view(request, project_uid):
    """Overview of a single project.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    context = {
        'project': project,
        'init_js_data': init_js_data,
        'tab_root': TAB_ROOT__DATA
    }
    return render(request, 'project.html', context)


@login_required
def project_delete(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    project.delete()
    response_data = {'redirect': '/'}
    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


@login_required
def tab_root_data(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA
    }
    return render(request, 'project.html', context)


VALID_ANALYZE_SUB_VIEWS = set([
    'variants',
    'sets',
    'genes'
])

@login_required
def tab_root_analyze(request, project_uid, ref_genome_uid=None, sub_view=None):
    context = {}

    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    # Initial javascript data.
    init_js_data = {
        'project': adapt_model_instance_to_frontend(project)
    }

    ref_genome = None
    if ref_genome_uid is not None:
        ref_genome = ReferenceGenome.objects.get(project=project,
                uid=ref_genome_uid)
        init_js_data['refGenome'] = adapt_model_instance_to_frontend(ref_genome)
        context['selected_ref_genome_uid'] = ref_genome_uid

    if sub_view is not None:
        if not sub_view in VALID_ANALYZE_SUB_VIEWS:
            raise Http404
        init_js_data['subView'] = sub_view
        context['sub_view'] = sub_view

    context.update({
        'project': project,
        'init_js_data': json.dumps(init_js_data),
        'tab_root': TAB_ROOT__ANALYZE,
    })

    return render(request, 'tab_root_analyze.html', context)


@login_required
def reference_genome_list_view(request, project_uid):
    """Shows the ReferenceGenomes and handles creating new
    ReferenceGenomes when requested.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    error_string = None

    # If a POST, then we are creating a new genome.
    if request.method == 'POST':
        # TODO: Add more inforative error handling.
        try:
            if request.POST['postType'] == 'modalFromFile':
                import_reference_genome_from_local_file(
                        project,
                        request.POST['refGenomeLabel'],
                        request.POST['refGenomeFileLocation'],
                        request.POST['importFileFormat'])
            elif request.POST['postType'] == 'modalFromNCBI':
                import_reference_genome_from_ncbi(
                        project,
                        request.POST['refGenomeLabel'],
                        request.POST['refGenomeAccession'],
                        request.POST['importFileFormat'])
            else:
                raise InputError('POST has invalid or missing postType')
        except Exception as e:
            error_string = 'Import error: ' + str(e)

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'ref_genome_list_json': adapt_model_to_frontend(ReferenceGenome,
                {'project':project}),
        'error_string': error_string
    }

    return render(request, 'reference_genome_list.html', context)


@login_required
def reference_genome_view(request, project_uid, ref_genome_uid):
    """Overview of a single project.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    reference_genome = ReferenceGenome.objects.get(project=project,
            uid=ref_genome_uid)
    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'reference_genome': reference_genome,
        'jbrowse_link': reference_genome.get_client_jbrowse_link()
    }
    return render(request, 'reference_genome.html', context)


@login_required
def sample_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    error_string = None

    # If a POST, then we are creating a new genome.
    if request.method == 'POST':
        # TODO: Add more informative error handling
        try:
            import_samples_from_targets_file(
                    project,
                    request.FILES['targetsFile'])
        except Exception as e:
            error_string = 'Import error: ' + str(e)

    # Query the db for the samples for this project.
    sample_list = ExperimentSample.objects.filter(project=project)

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'sample_list': sample_list,
        'error_string': error_string
    }
    return render(request, 'sample_list.html', context)


@login_required
def sample_list_targets_template(request):
    """Let the user download a blank sample targets template as a tab
    separated values file (.tsv) so they can fill it in and upload
    it back to the server.
    """
    context = {}
    return render(request, 'sample_list_targets_template.tsv', context,
            content_type='text/tab-separated-values')


@login_required
def variant_set_upload_template(request):
    """Let the user download a blank variant set template as a blank
    VCF file to be filled in.
    """
    context = {}
    return render(request, 'variant_set_upload_template.vcf', context,
            content_type='text/tab-separated-values')


@login_required
def alignment_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'alignment_list_json': adapt_model_to_frontend(AlignmentGroup,
                {'reference_genome__project':project})
    }
    return render(request, 'alignment_list.html', context)


@login_required
def alignment_view(request, project_uid, alignment_group_uid):
    """View of a single AlignmentGroup.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    alignment_group = AlignmentGroup.objects.get(
            reference_genome__project=project, uid=alignment_group_uid)

    if request.POST:
        run_snp_calling_pipeline(alignment_group)
        return HttpResponse('ok')

    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(alignment_group)
    })
    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'alignment_group': alignment_group,
        'experiment_sample_to_alignment_list_json': adapt_model_to_frontend(
                ExperimentSampleToAlignment,
                {'alignment_group': alignment_group}),
        'init_js_data': init_js_data
    }
    return render(request, 'alignment.html', context)


@login_required
def alignment_create_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    if request.POST:
        # Parse the data from the request body.
        request_data = json.loads(request.body)

        # Make sure the required keys are present.
        REQUIRED_KEYS = ['name', 'refGenomeUidList', 'sampleUidList']
        if not all(key in request_data for key in REQUIRED_KEYS):
            return HttpResponseBadRequest("Invalid request. Missing keys.")

        # Parse the data and look up the relevant model instances.
        alignment_group_name = request_data['name']
        if not len(alignment_group_name):
            return HttpResponseBadRequest("Name required.")

        ref_genome_list = ReferenceGenome.objects.filter(
                project=project,
                uid__in=request_data['refGenomeUidList'])
        assert len(ref_genome_list) == len(request_data['refGenomeUidList'])
        if not len(ref_genome_list) > 0:
            return HttpResponseBadRequest(
                    "At least one reference genome required.")

        sample_list = ExperimentSample.objects.filter(
                project=project,
                uid__in=request_data['sampleUidList'])
        assert len(sample_list) == len(request_data['sampleUidList'])
        if not len(sample_list) > 0:
            return HttpResponseBadRequest("At least one sample required.")

        # Kick off alignments.
        try:
            create_alignment_groups_and_start_alignments(alignment_group_name,
                    ref_genome_list, sample_list)

            # Success. Return a redirect response.
            response_data = {
                'redirect': reverse(
                        'genome_designer.main.views.alignment_list_view',
                        args=(project.uid,)),
            }
        except Exception as e:
            response_data = {
                'error': str(e)
            }

        return HttpResponse(json.dumps(response_data),
                content_type='application/json')

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'samples_list_json': adapt_model_to_frontend(ExperimentSample,
                {'project':project}),
        'ref_genomes_list_json': adapt_model_to_frontend(ReferenceGenome,
                {'project':project})
    }
    return render(request, 'alignment_create.html', context)

@login_required
def sample_alignment_error_view(request, project_uid, alignment_group_uid,
        sample_alignment_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    sample_alignment = ExperimentSampleToAlignment.objects.get(
            uid=sample_alignment_uid,
            alignment_group__reference_genome__project__uid=project_uid)

    # Get the path of the error file.
    data_dir = sample_alignment.get_model_data_dir()
    error_file_dir = os.path.join(data_dir, 'bwa_align.error')
    if os.path.exists(error_file_dir):
        with open(error_file_dir) as fh:
            raw_data = fh.read()
    else:
        raw_data = 'undefined'
    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'raw_data': raw_data
    }
    return render(request, 'sample_alignment_error_view.html', context)


@login_required
def variant_set_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    error_string = ''

    # If a POST, then we are creating a new variant set.
    if request.method == 'POST':
        # Common validation.
        # TODO: Most of this should be validated client-side.
        if (not 'refGenomeID' in request.POST or
                request.POST['refGenomeID'] == ''):
            error_string = 'No reference genome selected.'
        elif (not 'variantSetName' in request.POST or
                request.POST['variantSetName'] == ''):
            error_string = 'No variant set name.'

        if not error_string:
            ref_genome_uid = request.POST['refGenomeID']
            variant_set_name = request.POST['variantSetName']

            # Handling depending on which form was submitted.
            if request.POST['create-set-type'] == 'from-file':
                error_string = _create_variant_set_from_file(request, project,
                        ref_genome_uid, variant_set_name)
            elif request.POST['create-set-type'] == 'empty':
                error_string = _create_variant_set_empty(project,
                        ref_genome_uid, variant_set_name)
            else:
                return HttpResponseBadRequest("Invalid request.")

    # Grab all the ReferenceGenomes for this project
    # (for choosing which ref genome a new variant set goes into).
    ref_genome_list = ReferenceGenome.objects.filter(project=project)

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'ref_genome_list': ref_genome_list,
        'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                {'reference_genome__project':project}),
        'error_string': error_string
    }

    return render(request, 'variant_set_list.html', context)


def _create_variant_set_from_file(request, project, ref_genome_uid,
        variant_set_name):
    """Creates a variant set from file.

    Returns:
        A string indicating any errors that occurred. If no errors, then
        the empty string.
    """
    error_string = ''

    path = default_storage.save('tmp/tmp_varset.vcf',
            ContentFile(request.FILES['vcfFile'].read()))
    variant_set_file = os.path.join(settings.MEDIA_ROOT, path)

    try:
        ref_genome = ReferenceGenome.objects.get(project=project,
                uid=ref_genome_uid)
        import_variant_set_from_vcf(ref_genome, variant_set_name,
                variant_set_file)
    except Exception as e:
        error_string = 'Import error: ' + str(e)
    finally:
        os.remove(variant_set_file)

    return error_string


def _create_variant_set_empty(project, ref_genome_uid, variant_set_name):
    """Creates an empty variant set.

    Returns:
        A string indicating any errors that occurred. If no errors, then
        the empty string.
    """
    error_string = ''

    ref_genome = ReferenceGenome.objects.get(
            project=project,
            uid=ref_genome_uid)

    variant_set = VariantSet.objects.create(
            reference_genome=ref_genome,
            label=variant_set_name)

    return error_string


@login_required
def variant_set_view(request, project_uid, variant_set_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    variant_set = VariantSet.objects.get(
            reference_genome__project=project,
            uid=variant_set_uid)

    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(variant_set)
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'variant_set': variant_set,
        'variant_to_variant_set_json': adapt_model_to_frontend(
                VariantToVariantSet,
                {'variant_set': variant_set}),
        'init_js_data': init_js_data
    }

    return render(request, 'variant_set.html', context)


@login_required
def single_variant_view(request, project_uid, ref_genome_uid, variant_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    reference_genome = ReferenceGenome.objects.get(project=project,
        uid=ref_genome_uid)

    variant = get_object_or_404(Variant,
            uid=variant_uid,
            reference_genome=reference_genome)

    melted_variant_list = MeltedVariantView.variant_as_melted_list(variant)

    fe_melted_variant_list = adapt_model_or_modelview_list_to_frontend(
            melted_variant_list,
            variant_key_map=reference_genome.variant_key_map)

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'variant': variant,
        'melted_variant_list': fe_melted_variant_list
    }

    return render(request, 'single_variant_view.html', context)
