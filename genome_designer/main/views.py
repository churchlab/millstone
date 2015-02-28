"""
Views of pages.

NOTE: Put new Ajax-only actions into main/xhr_handlers.py.
"""

import json
import os
import re
import urllib

from django.conf import settings
from django.contrib.auth import login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.urlresolvers import reverse
from django.http import Http404
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.shortcuts import render
from registration.backends.simple.views import RegistrationView

from main.adapters import adapt_model_instance_to_frontend
from main.adapters import adapt_model_to_frontend
from main.models import AlignmentGroup
from main.models import Dataset
from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import VariantSet
from main.model_utils import get_dataset_with_type
from pipeline.pipeline_runner import run_pipeline
from utils.import_util import import_variant_set_from_vcf
from utils.jbrowse_util import compile_tracklist_json


# Tags used to indicate which tab we are on.
TAB_ROOT__DATA = 'DATA'
TAB_ROOT__ANALYZE = 'ANALYZE'


def home_view(request):
    """The main landing page.
    """
    context = {}
    return render(request, 'home.html', context)


def demo_splash_view(request):
    """The main landing page.
    """
    context = {}
    return render(request, 'demo_splash.html', context)


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

    context = {}

    if request.POST:
        try:
            project_name = request.POST.get('title')

            existing_proj_count = Project.objects.filter(
                    owner=request.user.get_profile(),
                    title=project_name).count()

            assert existing_proj_count == 0, (
                'Project with that name already exists.')

            project = Project.objects.create(
                    owner=request.user.get_profile(),
                    title=project_name)

            # Redirect to "Create Alignment" page.
            redirect_url = '%s?is_new_project=1' % reverse(
                    'main.views.alignment_create_view',
                    args=(project.uid,))
            return HttpResponseRedirect(redirect_url)

        except Exception as e:
            context['error_string'] = str(e)

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


VALID_ANALYZE_SUB_VIEWS = set([
    'variants',
    'sets',
    'genes'
])


@login_required
def tab_root_analyze(request, project_uid, alignment_group_uid=None, sub_view=None):
    """The analyze view. Subviews are loaded in Javascript.
    """
    context = {}

    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    # Initial javascript data to establish initial state.
    init_js_data = {
        'project': adapt_model_instance_to_frontend(project)
    }

    alignment_group = None
    if alignment_group_uid is not None:
        alignment_group = get_object_or_404(AlignmentGroup,
                uid=alignment_group_uid,
                reference_genome__project__owner=request.user.get_profile())
        init_js_data['alignmentGroup'] = adapt_model_instance_to_frontend(
                alignment_group)
        context['active_alignment_group'] = alignment_group

        ref_genome = alignment_group.reference_genome
        init_js_data['refGenome'] = adapt_model_instance_to_frontend(
                ref_genome)
        context['active_ref_genome'] = ref_genome

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
def genome_finish_view(request, project_uid):
    """Displays reference genomes on the opening page of
    the genome finish tab
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'init_js_data': init_js_data,
    }

    return render(request, 'genome_finish.html', context)

@login_required
def reference_genome_list_view(request, project_uid):
    """Shows the ReferenceGenomes and handles creating new
    ReferenceGenomes when requested.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'init_js_data': init_js_data,
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

    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'init_js_data': init_js_data,
        'error_string': error_string
    }
    return render(request, 'sample_list.html', context)

@login_required
def fastqc_view(request, project_uid, sample_uid, read_num):

    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    sample = get_object_or_404(ExperimentSample, project=project,
            uid=sample_uid)

    if int(read_num) == 1:
        dataset_type = Dataset.TYPE.FASTQC1_HTML
    elif int(read_num) == 2:
        dataset_type = Dataset.TYPE.FASTQC2_HTML
    else:
        raise Exception('Read number must be 1 or 2')

    fastqc_dataset = get_dataset_with_type(sample, dataset_type)

    response = HttpResponse(mimetype="text/html")
    for line in open(fastqc_dataset.get_absolute_location()):
        response.write(line)
    return response


@login_required
def alignment_list_view(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'init_js_data': init_js_data,
    }

    return render(request, 'alignment_list.html', context)


@login_required
def alignment_view(request, project_uid, alignment_group_uid):
    """View of a single AlignmentGroup.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    alignment_group = get_object_or_404(AlignmentGroup,
            reference_genome__project=project, uid=alignment_group_uid)

    if request.POST:
        return _rerun_alignment(alignment_group)

    # Initial javascript data.
    init_js_data = json.dumps({
        'project': adapt_model_instance_to_frontend(project),
        'alignment_group': adapt_model_instance_to_frontend(alignment_group)
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


def _rerun_alignment(alignment_group):
    """Re-runs existing alignment.
    """
    run_pipeline(alignment_group.label, alignment_group.reference_genome,
            alignment_group.get_samples())
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
def alignment_create_view(request, project_uid):
    """Displays the view for creating a new alignment.

    Also handles POST request to actually creating the alignment.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)

    # Validate and kick off new alignment.
    if request.POST:
        return _start_new_alignment(request, project)

    # Show the data that renders the create view.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    is_new_project = bool(request.GET.get('is_new_project', False))

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'init_js_data': init_js_data,
        'is_new_project': is_new_project
    }
    return render(request, 'alignment_create.html', context)


def _start_new_alignment(request, project):
    """Delegate function that handles logic of kicking off alignment.
    """
    # Parse the data from the request body.
    request_data = json.loads(request.body)

    # Make sure the required keys are present.
    REQUIRED_KEYS = [
            'name', 'refGenomeUidList', 'sampleUidList', 'skipHetOnly',
            'callAsHaploid']

    if not all(key in request_data for key in REQUIRED_KEYS):
        return HttpResponseBadRequest("Invalid request. Missing keys.")

    try:
        # Parse the data and look up the relevant model instances.
        alignment_group_name = request_data['name']
        assert len(alignment_group_name), "Name required."

        ref_genome_list = ReferenceGenome.objects.filter(
                project=project,
                uid__in=request_data['refGenomeUidList'])
        assert (len(ref_genome_list) ==
                len(request_data['refGenomeUidList'])), (
                        "Invalid reference genome uid(s).")
        assert len(ref_genome_list) == 1, (
            "Exactly one reference genome must be provided.")
        ref_genome = ref_genome_list[0]

        # Make sure AlignmentGroup has a unique name, because run_pipeline
        # will re-use an alignment based on label, reference genome,
        # aligner. We are currently hard-coding the aligner to BWA.
        assert AlignmentGroup.objects.filter(
                label=alignment_group_name,
                reference_genome=ref_genome).count() == 0, (
                        "Please pick unique alignment name.")

        sample_list = ExperimentSample.objects.filter(
                project=project,
                uid__in=request_data['sampleUidList'])
        assert len(sample_list) == len(request_data['sampleUidList']), (
                "Invalid expeirment sample uid(s).")
        assert len(sample_list) > 0, "At least one sample required."

        # Populate alignment options.
        alignment_options = dict()
        if request_data['skipHetOnly']:
            alignment_options['skip_het_only'] = True

        if request_data['callAsHaploid']:
            alignment_options['call_as_haploid'] = True

        # Kick off alignments.
        run_pipeline(
                alignment_group_name,
                ref_genome, sample_list,
                alignment_options=alignment_options)

        # Success. Return a redirect response.
        response_data = {
            'redirect': reverse(
                    'main.views.alignment_list_view',
                    args=(project.uid,)),
        }
    except Exception as e:
        response_data = {
            'error': str(e)
        }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


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

    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': adapt_model_instance_to_frontend(project)
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'ref_genome_list' : ref_genome_list,
        'init_js_data' : init_js_data,
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
    """Data view for a single VariantSet which includes options for taking
    actions on the VariantSet.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    variant_set = get_object_or_404(VariantSet,
            reference_genome__project=project,
            uid=variant_set_uid)

    # Initial javascript data.
    init_js_data = json.dumps({
        'entity': {
            'uid': variant_set.uid,
            'project': adapt_model_instance_to_frontend(project),
            'refGenomeUid': variant_set.reference_genome.uid,
        }
    })

    context = {
        'project': project,
        'tab_root': TAB_ROOT__DATA,
        'variant_set': variant_set,
        'init_js_data': init_js_data,
        'is_print_mage_oligos_enabled':
                settings.FLAG__PRINT_MAGE_OLIGOS_ENABLED,
        'is_generate_new_reference_genome_enabled':
                settings.FLAG__GENERATE_NEW_REFERENCE_GENOME_ENABLED
    }

    return render(request, 'variant_set.html', context)


@login_required
def compile_jbrowse_and_redirect(request):
    """Compiles the jbrowse tracklist and redirects to jbrowse
    """
    # First, grab the data string and get the project and ref genome from it
    data_string = request.GET.get('data')
    regexp_str = (r'/jbrowse/gd_data/projects/(?P<project_uid>\w+)' +
            r'/ref_genomes/(?P<ref_genome_uid>\w+)/jbrowse')
    uid_match = re.match(regexp_str, data_string)
    assert uid_match, "Incorrect URL passed in data key to redirect_jbrowse"

    reference_genome = get_object_or_404(ReferenceGenome,
        project__owner=request.user.get_profile(),
        uid=uid_match.group('ref_genome_uid'))

    # Recompile the tracklist from components and symlink subdirs
    compile_tracklist_json(reference_genome)

    # Finally, pass the GET along to the jbrowse static page.
    get_values = urllib.urlencode(request.GET).replace('+', '%20')

    maybe_force_nginx = ''
    if settings.DEBUG_FORCE_JBROWSE_NGINX:
        maybe_force_nginx = 'http://localhost'
    full_redirect_url = (maybe_force_nginx + '/jbrowse/index.html' + '?' +
            get_values)

    return HttpResponseRedirect(full_redirect_url)


class RegistrationViewWrapper(RegistrationView):
    """
    For now, if there are no users present, allow registration.
    Once the first user is created, disallow registration.
    """

    def registration_allowed(self, request):

        if not User.objects.count():
            return True
        return super(RegistrationViewWrapper, self).registration_allowed(
                request)

    def get_success_url(self, request, user):
        '''
        Log in the new user and take them to the project list.
        '''
        login(request, user)
        return '/'


if settings.RUNNING_ON_EC2:
    from boto.utils import get_instance_metadata
    @login_required
    def ec2_info_view(request):
        """
        boto.utils.get_instance_metadata() will block if not running on EC2.
        """
        m = get_instance_metadata()

        password_path = os.path.join(settings.PWD, "../password.txt")
        if os.path.isfile(password_path):
            with open(password_path) as f:
                password = f.read().strip()
        else:
            password = "Not found from '%s'" % password_path

        context = {
            'hostname': m['public-hostname'],
            'instance_type': m['instance-type'],
            'instance_id': m['instance-id'],
            'password': password
        }

        return render(request, 'ec2_info_view.html', context)
