import json

from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.shortcuts import render

from main.adapters import adapt_model_instance_to_frontend
from main.adapters import adapt_model_to_frontend
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
from scripts.import_util import import_samples_from_targets_file
from scripts.snp_callers import run_snp_calling_pipeline

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
def project_view(request, project_uid):
    """Overview of a single project.
    """
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project,
    }
    return render(request, 'project.html', context)


@login_required
def reference_genome_list_view(request, project_uid):
    """Shows the ReferenceGenomes and handles creating new
    ReferenceGenomes when requested.
    """
    project = Project.objects.get(uid=project_uid)
    error_string = None

    # If a POST, then we are creating a new genome.
    if request.method == 'POST':
        # TODO: Add more inforative error handling.
        try:
            import_reference_genome_from_local_file(
                    project,
                    request.POST['refGenomeLabel'],
                    request.POST['refGenomeFileLocation'],
                    request.POST['importFileFormat'])
        except Exception as e:
            error_string = 'Import error: ' + str(e)

    # Grab all the ReferenceGenomes for this project.
    ref_genome_list = ReferenceGenome.objects.filter(project=project)

    # Adapt the backend objects to the frontend format.
    fe_ref_genome_list = [{
        'label': obj.label,
        'num_chromosomes': obj.num_chromosomes,
        'total_size': obj.num_bases,
        'annotated': False,
        'parents': [],
        'children': [],
        'href': obj.get_href()
    } for obj in ref_genome_list]

    # HACK: Add some fake data for now.
    fe_ref_genome_list.append({
        'label': 'MG1655-fake',
        'num_chromosomes': '1',
        'total_size': '3.6 Mbp',
        'annotated': False,
        'parents': [],
        'children': ['C321D-fake', 'another'],
    })

    fe_ref_genome_list.append({
        'label': 'C321D-fake',
        'num_chromosomes': '1',
        'total_size': '3.6 Mbp',
        'annotated': True,
        'parents': ['MG1655-fake'],
        'children': [],
    })

    # (Re-)Render the page.
    context = {
        'project': project,
        'ref_genome_list': fe_ref_genome_list,
        'error_string': error_string
    }
    return render(request, 'reference_genome_list.html', context)


@login_required
def reference_genome_view(request, project_uid, ref_genome_uid):
    """Overview of a single project.
    """
    project = Project.objects.get(uid=project_uid)
    reference_genome = ReferenceGenome.objects.get(uid=ref_genome_uid)
    context = {
        'project': project,
        'reference_genome': reference_genome,
        'jbrowse_link': reference_genome.get_client_jbrowse_link()
    }
    return render(request, 'reference_genome.html', context)


@login_required
def sample_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    error_string = None

    # If a POST, then we are creating a new genome.
    if request.method == 'POST':
        # TODO: Add more inforative error handling
        print "Recieved POST!"
        try:
            import_samples_from_targets_file(
                    project,
                    request.FILES['targetsFile'])
        except Exception as e:
            error_string = 'Import error: ' + str(e)

    context = {
        'project': project,
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
def alignment_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)

    context = {
        'project': project,
        'alignment_list_json': adapt_model_to_frontend(AlignmentGroup,
                {'reference_genome__project':project})
    }
    return render(request, 'alignment_list.html', context)


@login_required
def alignment_view(request, project_uid, alignment_group_uid):
    """View of a single AlignmentGroup.
    """
    project = Project.objects.get(uid=project_uid)
    alignment_group = AlignmentGroup.objects.get(uid=alignment_group_uid)

    if request.POST:
        run_snp_calling_pipeline(alignment_group)
        return HttpResponse('ok')

    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(alignment_group)
    })
    context = {
        'project': project,
        'alignment_group': alignment_group,
        'experiment_sample_to_alignment_list_json': adapt_model_to_frontend(
                ExperimentSampleToAlignment,
                {'alignment_group': alignment_group}),
        'init_js_data': init_js_data
    }
    return render(request, 'alignment.html', context)


@login_required
def alignment_create_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)

    if request.POST:
        # Parse the data from the request body.
        request_data = json.loads(request.body)

        # Make sure the required keys are present.
        REQUIRED_KEYS = ['refGenomeUidList', 'sampleUidList']
        if not all(key in request_data for key in REQUIRED_KEYS):
            return HttpResponseBadRequest("Invalid request. Missing keys.")

        # Parse the data and look up the relevant model instances.

        ref_genome_list = ReferenceGenome.objects.filter(
                uid__in=request_data['refGenomeUidList'])
        assert len(ref_genome_list) == len(request_data['refGenomeUidList'])
        if not len(ref_genome_list) > 0:
            return HttpResponseBadRequest(
                    "At least one reference genome required.")

        sample_list = ExperimentSample.objects.filter(
                uid__in=request_data['sampleUidList'])
        assert len(sample_list) == len(request_data['sampleUidList'])
        if not len(sample_list) > 0:
            return HttpResponseBadRequest(
                    "At least sample required.")

        # Kick off alignments.
        # NOTE: Hard-coded test_models_only=True for now.
        create_alignment_groups_and_start_alignments(ref_genome_list,
                sample_list, test_models_only=False)

        # Success. Return a redirect response.
        response_data = {
            'redirect': reverse(
                    'genome_designer.main.views.alignment_list_view',
                    args=(project.uid,)),
        }
        return HttpResponse(json.dumps(response_data),
                content_type='application/json')

    context = {
        'project': project,
        'samples_list_json': adapt_model_to_frontend(ExperimentSample,
                {'project':project}),
        'ref_genomes_list_json': adapt_model_to_frontend(ReferenceGenome,
                {'project':project})
    }
    return render(request, 'alignment_create.html', context)


@login_required
def variant_set_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
        
    context = {
        'project': project,
        'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                {'reference_genome__project':project})
    }
    
    return render(request, 'variant_set_list.html', context)

@login_required
def variant_set_view(request, project_uid, variant_set_uid):
    project = Project.objects.get(uid=project_uid)
    variant_set = VariantSet.objects.get(uid=variant_set_uid)
    
    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(variant_set)
    })
    
    context = {
        'project': project,
        'variant_set': variant_set,
        'variant_to_variant_set_json': adapt_model_to_frontend(
                VariantToVariantSet,
                {'variant_set': variant_set}),
        'init_js_data': init_js_data
    }
    
    return render(request, 'variant_set.html', context)


@login_required
def variant_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)

    # Fetch the list of variants and render it into the dom as json.
    # The data will be displayed to the user via the javascript DataTables
    # component.
    context = {
       'project': project,
       'variant_list_json': adapt_model_to_frontend(Variant,
            {'reference_genome__project':project})
    }

    return render(request, 'variant_list.html', context)


@login_required
def gene_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project
    }
    return render(request, 'gene_list.html', context)


@login_required
def goterm_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project
    }
    return render(request, 'goterm_list.html', context)
