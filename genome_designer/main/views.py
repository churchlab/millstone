import json

from django.contrib.auth.decorators import login_required
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.urlresolvers import reverse
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.shortcuts import render
import os

from main.adapters import adapt_model_instance_to_frontend
from main.adapters import adapt_model_to_frontend
from main.forms import ProjectForm
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
from scripts.import_util import import_variant_set_from_vcf
from scripts.snp_callers import run_snp_calling_pipeline
from scripts.variant_sets import add_or_remove_variants_from_set
import settings

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

    context = {
        'project': project,
    }
    return render(request, 'project.html', context)


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
            import_reference_genome_from_local_file(
                    project,
                    request.POST['refGenomeLabel'],
                    request.POST['refGenomeFileLocation'],
                    request.POST['importFileFormat'])
        except Exception as e:
            error_string = 'Import error: ' + str(e)

    context = {
        'project': project,
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
        REQUIRED_KEYS = ['refGenomeUidList', 'sampleUidList']
        if not all(key in request_data for key in REQUIRED_KEYS):
            return HttpResponseBadRequest("Invalid request. Missing keys.")

        # Parse the data and look up the relevant model instances.

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
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    error_string = ''

    # If a POST, then we are creating a new variant set.
    if request.method == 'POST':
        # TODO: Add more informative error handling.

        # Save vcf file to disk temporarily first.
        path = default_storage.save('tmp/tmp_varset.vcf',
                ContentFile(request.FILES['vcfFile'].read()))
        variant_set_file = os.path.join(settings.MEDIA_ROOT, path)

        try:
            # First validate the request.
            # TODO: Most of this should be validated client-side.
            if not 'refGenomeID' in request.POST:
                error_string = 'No reference genome selected.'
            elif (not 'variantSetName' in request.POST or
                    request.POST['variantSetName'] == ''):
                error_string = 'No variant set name.'

            # If no error here, then continue.
            if not error_string:
                import_variant_set_from_vcf(
                        project,
                        request.POST['refGenomeID'],
                        request.POST['variantSetName'],
                        variant_set_file)
        except Exception as e:
            error_string = 'Import error: ' + str(e)
        finally:
            os.remove(variant_set_file)

    # Grab all the ReferenceGenomes for this project
    # (for choosing which ref genome a new variant set goes into).
    ref_genome_list = ReferenceGenome.objects.filter(project=project)

    # We only need the label and ID for every reference genome
    fe_ref_genome_list = [{
        'label': obj.label,
        'id': obj.id} for obj in ref_genome_list]

    context = {
        'project': project,
        'ref_genome_list': fe_ref_genome_list,
        'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                {'reference_genome__project':project}),
        'error_string': error_string
    }

    return render(request, 'variant_set_list.html', context)

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
        'variant_set': variant_set,
        'variant_to_variant_set_json': adapt_model_to_frontend(
                VariantToVariantSet,
                {'variant_set': variant_set}),
        'init_js_data': init_js_data
    }

    return render(request, 'variant_set.html', context)


@login_required
def variant_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    # The json data required to populate datables and dropdowns.
    # Can be called either by render or post response.
    def get_json_data():
        return {
           'variant_list_json': adapt_model_to_frontend(Variant,
                {'reference_genome__project':project}),
           'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                {'reference_genome__project':project})}

    context = {}

    # If we are returning data after a post and not redirecting.
    if request.POST:
        # Parse the data from the request body.
        request_data = json.loads(request.body)

        # Make sure the required keys are present.
        REQUIRED_KEYS = [
                'variantUidList',
                'variantSetAction',
                'variantSetUid']

        if not all(key in request_data for key in REQUIRED_KEYS):
            return HttpResponseBadRequest("Invalid request. Missing keys.")

        # Add or remove the variants to the set, as variantSetAction.
        context.update(add_or_remove_variants_from_set(**request_data))

        # Get new updated json info for datatables and dropdown.
        context.update(get_json_data())

        return HttpResponse(
                json.dumps(context),
                content_type='application/json')

    # If we are rendering a new view.
    else:
        context.update({'project': project})
        context.update(get_json_data())
        return render(request, 'variant_list.html', context)


@login_required
def gene_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    context = {
        'project': project
    }
    return render(request, 'gene_list.html', context)


@login_required
def goterm_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    context = {
        'project': project
    }
    return render(request, 'goterm_list.html', context)
