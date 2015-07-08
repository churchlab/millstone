"""
Methods that handle Ajax requests from the frontend.

This module was created in response to views.py getting quite big, and so a
reasonable separation point is to separate page actions from Ajax actions.
"""

import copy
import json
import os
from StringIO import StringIO
import tempfile

from Bio import SeqIO
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.servers.basehttp import FileWrapper
from django.core.urlresolvers import reverse
from django.http import Http404
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import StreamingHttpResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_GET
from django.views.decorators.http import require_POST

#from debug.debug_util import FakeException
from main.adapters import adapt_model_to_frontend
from main.adapters import adapt_experiment_samples_to_frontend
from main.exceptions import ValidationException
from main.model_utils import get_dataset_with_type
from main.model_views import adapt_gene_list_to_frontend
from main.model_views import get_all_fields
from main.model_views import adapt_variant_to_frontend
from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import Contig
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import SavedVariantFilterQuery
from main.models import VariantCallerCommonData
from main.models import VariantAlternate
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import S3File
from genome_finish import assembly
from genome_finish.insertion_placement import find_contig_insertion_site
from genome_finish.insertion_placement import place_cassette
from utils.combine_reference_genomes import combine_list_allformats
from utils.data_export_util import export_melted_variant_view
from utils.import_util import create_samples_from_row_data
from utils.import_util import create_sample_models_for_eventual_upload
from utils.import_util import import_reference_genome_from_local_file
from utils.import_util import import_reference_genome_from_ncbi
from utils.import_util import import_samples_from_targets_file
from utils.import_util import import_variant_set_from_vcf
from utils.optmage_util import ReplicationOriginParams
from utils.optmage_util import print_mage_oligos
from utils.reference_genome_maker_util import generate_new_reference_genome
from variants.common import determine_visible_field_names
from variants.filter_key_map_constants import MAP_KEY__ALTERNATE
from variants.filter_key_map_constants import MAP_KEY__COMMON_DATA
from variants.filter_key_map_constants import MAP_KEY__EVIDENCE
from variants.gene_query import lookup_genes
from variants.materialized_variant_filter import lookup_variants
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager
from variants.variant_sets import update_variant_in_set_memberships
from variants.variant_sets import update_variant_in_set_memberships__all_matching_filter

if settings.S3_ENABLED:
    from utils.import_util import parse_targets_file, import_reference_genome_from_s3, import_samples_from_s3
    from s3 import s3_get_string


@login_required
def project_delete(request, project_uid):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    project.delete()
    response_data = {'redirect': '/'}
    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


@login_required
@require_POST
def create_ref_genome_from_browser_upload(request):
    """Handle request to create ReferenceGenome from local file.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST['projectUid'])

    uploaded_file = request.FILES['refGenomeFile']

    # Save uploaded ReferenceGenome to temp file, passing the original filename
    # as the suffix for easier debug.
    if not os.path.exists(settings.TEMP_FILE_ROOT):
        os.mkdir(settings.TEMP_FILE_ROOT)
    _, temp_file_location = tempfile.mkstemp(
        suffix='_' + uploaded_file.name,
        dir=settings.TEMP_FILE_ROOT)

    with open(temp_file_location, 'w') as temp_file_fh:
        temp_file_fh.write(request.FILES['refGenomeFile'].read())

    error_string = ''
    try:
        import_reference_genome_from_local_file(
                project,
                request.POST['refGenomeLabel'],
                temp_file_location,
                request.POST['importFileFormat'],
                move=True)
    except Exception as e:
        error_string = str(e)

    result = {
        'error': error_string,
    }

    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def create_ref_genome_from_server_location(request):
    """Handle request to create ReferenceGenome from local file.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST['projectUid'])

    error_string = ''
    try:
        import_reference_genome_from_local_file(
                project,
                request.POST['refGenomeLabel'],
                request.POST['refGenomeFileLocation'],
                request.POST['importFileFormat'])
    except Exception as e:
        error_string = str(e)

    result = {
        'error': error_string,
    }

    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def create_ref_genome_from_ncbi(request):
    """Handle request to create ReferenceGenome from local file.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST['projectUid'])

    error_string = ''
    try:
        import_reference_genome_from_ncbi(
                project,
                request.POST['refGenomeLabel'],
                request.POST['refGenomeAccession'],
                request.POST['importFileFormat'])
    except Exception as e:
        error_string = str(e)

    result = {
        'error': error_string,
    }

    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def ref_genomes_delete(request):
    """Deletes ReferenceGenomes.
    """
    request_data = json.loads(request.body)
    ref_genome_uid_list = request_data.get('refGenomeUidList', [])
    if len(ref_genome_uid_list) == 0:
        raise Http404

    # First make sure all the samples belong to this user.
    ref_genomes_to_delete = ReferenceGenome.objects.filter(
            project__owner=request.user.get_profile(),
            uid__in=ref_genome_uid_list)
    if not len(ref_genomes_to_delete) == len(ref_genome_uid_list):
        raise Http404

    # Validation successful, delete.
    ref_genomes_to_delete.delete()

    # Return success response.
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_POST
def ref_genomes_concatenate(request):
    """Concatenates ReferenceGenomes.
    """
    request_data = json.loads(request.POST['data'])
    ref_genome_uid_list = request_data['refGenomeUidList']
    if len(ref_genome_uid_list) == 0:
        raise Http404
    new_genome_label = request_data['newGenomeLabel']
    if len(new_genome_label) == 0:
        raise Http404

    # First make sure all the samples belong to this user.
    ref_genomes_to_concatenate = ReferenceGenome.objects.filter(
            project__owner=request.user.get_profile(),
            uid__in=ref_genome_uid_list)
    if not len(ref_genomes_to_concatenate) == len(ref_genome_uid_list):
        raise Http404

    # Validation successful, concatenate.
    project = ref_genomes_to_concatenate[0].project
    combine_list_allformats(
            ref_genomes_to_concatenate, new_genome_label, project)

    # Return success response.
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_GET
def ref_genomes_download(request):
    """Downloads requested fasta/genbank file
    """
    file_format = request.GET['file_format']
    reference_genome = get_object_or_404(
            ReferenceGenome, uid=request.GET['reference_genome_uid'])
    if file_format == 'fasta':
        file_path = reference_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
        file_name = '.'.join([reference_genome.label, 'fa'])
    elif file_format == 'genbank':
        file_path = reference_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()
        file_name = '.'.join([reference_genome.label, 'gb'])
    else:
        raise Http404

    wrapper = FileWrapper(file(file_path))
    response = StreamingHttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="{0}"'.format(
            file_name)
    response['Content-Length'] = os.path.getsize(file_path)

    return response


@login_required
@require_GET
def contigs_download(request):
    """Downloads fasta file of contig sequence
    """
    contig = get_object_or_404(
            Contig, uid=request.GET['contig_uid'])

    file_path = contig.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    file_name = '.'.join([contig.label, 'fa'])

    wrapper = FileWrapper(file(file_path))
    response = StreamingHttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="{0}"'.format(
            file_name)
    response['Content-Length'] = os.path.getsize(file_path)

    return response


@login_required
@require_POST
def contigs_delete(request):
    """Deletes ReferenceGenomes.
    """
    request_data = json.loads(request.body)
    contig_uid_list = request_data.get('contigUidList', [])
    if len(contig_uid_list) == 0:
        raise Http404

    # First make sure all the samples belong to this user.
    contigs_to_delete = Contig.objects.filter(
            parent_reference_genome__project__owner=(
                    request.user.get_profile()),
            uid__in=contig_uid_list)
    if not len(contigs_to_delete) == len(contig_uid_list):
        raise Http404

    # Validation successful, delete.
    contigs_to_delete.delete()

    # Return success response.
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_POST
def variant_sets_delete(request):
    """Deletes a list of variant sets.
    """
    request_data = json.loads(request.body)
    variant_set_uid_list = request_data.get('variantSetUidList')

    # First make sure all the sets belong to this user.
    variant_sets_to_delete = VariantSet.objects.filter(
        reference_genome__project__owner=request.user.get_profile(),
        uid__in=variant_set_uid_list)

    if not len(variant_sets_to_delete) == len(variant_set_uid_list):
        raise Http404

    # Validation succcessful, delete
    variant_sets_to_delete.delete()

    # Return success response
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_POST
def save_variant_filter(request):
    # Read params / validate.
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST['projectUid'])
    filter_text = request.POST.get('filterText', '')
    if not filter_text:
        raise Http404("Nothing to save.")

    # Get or create the new filter.
    svfq, _ = SavedVariantFilterQuery.objects.get_or_create(
            owner=project.owner,
            text=filter_text)

    # Return new element in response so it can be rendered.
    return HttpResponse(json.dumps({
        'savedFilter': {
            'uid': svfq.uid,
            'text': svfq.text
        }
    }), content_type='application/json')


@login_required
@require_POST
def delete_variant_filter(request):
    # Read params / validate.
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST['projectUid'])
    filter_uid = request.POST.get('uid')
    if not filter_uid:
        return HttpResponseBadRequest();

    # Get or create the new filter.
    SavedVariantFilterQuery.objects.get(
            owner=project.owner, uid=filter_uid).delete()

    # Return new element in response so it can be rendered.
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_POST
def upload_single_sample(request):
    # Read params / validate.
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST['projectUid'])
    if not 'fastq1' in request.FILES:
        raise Http404
    sample_label = request.POST.get('sampleLabel', None)
    if not sample_label:
        raise Http404

    # Save uploaded Samples to temp location.
    if not os.path.exists(settings.TEMP_FILE_ROOT):
        os.mkdir(settings.TEMP_FILE_ROOT)
    fastq1_uploaded_file = request.FILES['fastq1']
    if not os.path.exists(settings.TEMP_FILE_ROOT):
        os.mkdir(settings.TEMP_FILE_ROOT)
    _, fastq1_temp_file_location = tempfile.mkstemp(
        suffix='_' + fastq1_uploaded_file.name,
        dir=settings.TEMP_FILE_ROOT)
    with open(fastq1_temp_file_location, 'w') as temp_file_fh:
        temp_file_fh.write(fastq1_uploaded_file.read())

    # Maybe handle fastq2.
    fastq2_uploaded_file = None
    if 'fastq2' in request.FILES:
        fastq2_uploaded_file = request.FILES['fastq2']
        _, fastq2_temp_file_location = tempfile.mkstemp(
            suffix='_' + fastq2_uploaded_file.name,
            dir=settings.TEMP_FILE_ROOT)
        with open(fastq2_temp_file_location, 'w') as temp_file_fh:
            temp_file_fh.write(fastq2_uploaded_file.read())

    result = {}

    # Create the data structure that the util expects and create samples.
    try:
        data_source_list = [{
            'Sample_Name': sample_label,
            'Read_1_Path': fastq1_temp_file_location,
            'Read_2_Path': fastq2_temp_file_location
        }]
        create_samples_from_row_data(project, data_source_list, move=False)
    except Exception as e:
        result['error'] = str(e)

    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def create_samples_from_server_location(request):
    """Handle request to create ReferenceGenome from local file.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST.get('projectUid', ''))

    try:
        import_samples_from_targets_file(
                project,
                request.FILES['targetsFile'])
    except Exception as e:
        result = {
            'error': str(e)
        }
        return HttpResponse(json.dumps(result),
                content_type='application/json')

    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_GET
def get_samples_awaiting_upload(request):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.GET.get('projectUid', ''))
    existing_sample_dataset_filename_list = [
            os.path.split(ds.filesystem_location)[1]
            for ds in Dataset.objects.filter(
                    experimentsample__project=project,
                    status=Dataset.STATUS.AWAITING_UPLOAD)]
    result = {
        'sampleFilenameList': existing_sample_dataset_filename_list
    }
    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def samples_upload_through_browser_template(request):
    """Handle request to create ReferenceGenome from local file.
    """
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST.get('projectUid', ''))

    try:
        template_file = request.FILES['file']
    except:
        result = {
            'error': 'Problem receiving file in request.'
        }
        return HttpResponse(json.dumps(result),
                content_type='application/json')

    try:
        create_sample_models_for_eventual_upload(project, template_file)
    except ValidationException as e:
        result = {
            'error': str(e)
        }
        return HttpResponse(json.dumps(result),
                content_type='application/json')

    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_POST
def samples_upload_through_browser_sample_data(request):
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=request.POST.get('projectUid', ''))

    # Grab the file from the request.
    uploaded_file = request.FILES['file']

    # Find the Dataset that matches the filename, validating at the same time.
    experiment_sample_datasets_in_project = Dataset.objects.filter(
            experimentsample__project=project)
    datasets_matching_project_and_filename = []
    for ds in experiment_sample_datasets_in_project:
        expected_filename = os.path.split(ds.filesystem_location)[1]
        if expected_filename == uploaded_file.name:
            datasets_matching_project_and_filename.append(ds)
    if len(datasets_matching_project_and_filename) == 0:
        result = {
            'error': 'UPLOAD ERROR: '
                     'Unexpected filename. Are you sure it\'s correct?'
        }
        return HttpResponse(json.dumps(result),
                content_type='application/json')

    # If this occurs, this is a bug. The upload should prevent Datasets with
    # the same filename for a particular project.
    assert len(datasets_matching_project_and_filename) == 1, (
            "Each Dataset must have a unique name.")

    # Identify the copy destination.
    dataset = datasets_matching_project_and_filename[0]
    copy_dest = dataset.get_absolute_location()

    # Copy the file in chunks.
    # TODO: Understand this better. Probably need error handling.
    with open(copy_dest, 'w') as dest_fh:
        for chunk in uploaded_file.chunks():
            dest_fh.write(chunk)

    # Update the status.
    dataset.status = Dataset.STATUS.READY
    dataset.save(update_fields=['status'])

    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_POST
def samples_delete(request):
    """Deletes ExperimentSamples that are not part of an AlignmentGroup.
    """
    request_data = json.loads(request.body)
    sample_uid_list = request_data.get('sampleUidList', [])
    if len(sample_uid_list) == 0:
        raise Http404

    # First make sure all the samples belong to this user.
    samples_to_delete = ExperimentSample.objects.filter(
            project__owner=request.user.get_profile(),
            uid__in=sample_uid_list)
    if not len(samples_to_delete) == len(sample_uid_list):
        raise Http404

    # Next, make sure none of these samples are part of an AlignmentGroup.
    samples_associated_with_alignment = []
    for sample in samples_to_delete:
        if (ExperimentSampleToAlignment.objects.filter(
                experiment_sample=sample).count() > 0):
            samples_associated_with_alignment.append(sample)
    if len(samples_associated_with_alignment) > 0:
        affected_samples = ', '.join([
                s.label for s in samples_associated_with_alignment])
        error_string = (
                '%s associated with an alignment. You must delete '
                'all related alignments for a sample before deleting it.' % (
                        affected_samples))
        result = {
            'error': error_string,
        }
        return HttpResponse(json.dumps(result), content_type='application/json')

    # Validation successful, delete.
    samples_to_delete.delete()

    # Return success response.
    return HttpResponse(json.dumps({}), content_type='application/json')


# Key in the GET params containing the string for filtering the variants.
VARIANT_LIST_REQUEST_KEY__FILTER_STRING = 'variantFilterString'
VARIANT_LIST_REQUEST_KEY__PROJECT_UID = 'projectUid'
VARIANT_LIST_REQUEST_KEY__REF_GENOME_UID = 'refGenomeUid'

VARIANT_LIST_RESPONSE_KEY__LIST = 'variant_list_json'
VARIANT_LIST_RESPONSE_KEY__TOTAL = 'num_total_variants'
VARIANT_LIST_RESPONSE_KEY__SET_LIST = 'variant_set_list_json'
VARIANT_LIST_RESPONSE_KEY__KEY_MAP = 'variant_key_filter_map_json'
VARIANT_LIST_RESPONSE_KEY__ERROR = 'error'

# Uncomment this and @profile statement to profile. This is the entry point
# to a monster SQL call so leaving this debugging code here commented out is
# useful.
#from debug.profiler import profile
#@profile('profile.log')
@login_required
@require_GET
def get_variant_list(request):
    """Returns a list of Variants, filtered by any filter parameters contained
    in the request.
    """
    # Parse the GET params.
    ref_genome_uid = request.GET.get('refGenomeUid')
    project_uid = request.GET.get('projectUid')
    maybe_alignment_group_uid = request.GET.get('alignmentGroupUid', None)

    # Get models and verify permissions.
    reference_genome = get_object_or_404(ReferenceGenome,
            project__uid=project_uid, uid=ref_genome_uid)
    if maybe_alignment_group_uid:
        alignment_group = get_object_or_404(AlignmentGroup,
                reference_genome=reference_genome,
                uid=maybe_alignment_group_uid)
    else:
        alignment_group = None

    # Dictionary to hold all query specific parameters.
    query_args = {}

    # Get inputs to perform the query for Variants data.
    # TODO: Combine with saved filter string.
    query_args['filter_string'] = request.GET.get(
            VARIANT_LIST_REQUEST_KEY__FILTER_STRING, '')

    # Determine whether melted or cast view.
    query_args['is_melted'] = request.GET.get('melt', 0) == '1'

    # Get optional column to sort by.
    # TODO shouldn't cast a client parameter to int outside of try-catch.
    query_args['sortCol'] = int(request.GET.get('iSortCol_0', 0))
    query_args['sort_by_direction'] = request.GET.get('sSortDir_0', 'asc')

    # Want all results listed, so set count_only to false.
    query_args['count_only'] = False

    # Pagination.
    query_args['pagination_start'] = int(request.GET.get('iDisplayStart', 0))
    query_args['pagination_len'] = int(request.GET.get('iDisplayLength', 100))

    # Any exception from here should be caused by a malformed query from the
    # user and the data should return an error string, rather than throw a 500.
    # Of course, it is possible that we have our bugs right now so devs should
    # be wary of this big try-except.
    try:
        field_select_keys = json.loads(request.GET.get(
                VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS, json.dumps([])))
        query_args['visible_key_names'] = determine_visible_field_names(
                field_select_keys, query_args['filter_string'],
                reference_genome)

        if query_args['sortCol']:  # 1 indexed; 0 means no sort column
            all_fields = get_all_fields(
                    reference_genome, query_args['visible_key_names'],
                    melted=query_args['is_melted'])
            # Get rid of hidden fields for sorting consideration.
            all_fields = [field for field in all_fields
                if not ('hide' in field and field['hide'])]
            if query_args['sortCol'] <= len(all_fields):
                query_args['sort_by_column'] = \
                    all_fields[query_args['sortCol'] - 1]['field']
            else:
                query_args['sort_by_column'] = ''
        else:
            query_args['sort_by_column'] = ''

        # Get the list of Variants (or melted representation) to display.
        lookup_variant_result = lookup_variants(query_args, reference_genome,
                alignment_group=alignment_group)
        variant_list = lookup_variant_result.result_list
        num_total_variants = lookup_variant_result.num_total_variants

        # Adapt the Variants to display for the frontend.
        variant_list_json = adapt_variant_to_frontend(variant_list,
                reference_genome, query_args['visible_key_names'],
                melted=query_args['is_melted'])
        
        # Get all VariantSets that exist for this ReferenceGenome.
        variant_set_list = VariantSet.objects.filter(
                reference_genome=reference_genome)

        # Query the keys valid for ReferenceGenome, and mark the ones that
        # will be displayed so that the checkmarks in the visible field select
        # are pre-filled in case the user wishes to change these.
        variant_key_map_with_active_fields_marked = copy.deepcopy(
                reference_genome.variant_key_map)
        _mark_active_keys_in_variant_key_map(
                variant_key_map_with_active_fields_marked,
                query_args['visible_key_names'])

        # Package up the response.
        response_data = {
            VARIANT_LIST_RESPONSE_KEY__LIST: variant_list_json,
            VARIANT_LIST_RESPONSE_KEY__TOTAL: num_total_variants,
            VARIANT_LIST_RESPONSE_KEY__SET_LIST: adapt_model_to_frontend(VariantSet,
                    obj_list=variant_set_list),
            VARIANT_LIST_RESPONSE_KEY__KEY_MAP: json.dumps(
                    variant_key_map_with_active_fields_marked)
        }
    # Toggle which of the following exceptions is commented for debugging.
    # except FakeException as e:
    except Exception as e:
        # TODO: More readable error reporting.
        exception_as_string = str(type(e)) + ' ' + str(e)
        response_data = {
            VARIANT_LIST_RESPONSE_KEY__ERROR: exception_as_string
        }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS = 'visibleKeyNames'


def _mark_active_keys_in_variant_key_map(variant_key_map, visible_key_names):
    """Mutates variant_key_map to mark fields that should be active based
    on model class defaults.
    """
    # In the current implementation, we mark the fields that are included
    # in the relevant models' get_field_order() method.

    def _update_model_class_key_map(model_class, variant_key_submap):
        """Helper method."""
        for key in visible_key_names:
            if key in variant_key_submap:
                variant_key_submap[key]['checked'] = True

        # TODO: Do we want to bring back this old default?
        # default_keys = [el['field'] for el in model_class.get_field_order()]
        # for key in default_keys:
        #     if key in variant_key_submap:
        #         variant_key_submap[key]['checked'] = True

    _update_model_class_key_map(VariantCallerCommonData,
            variant_key_map[MAP_KEY__COMMON_DATA])

    _update_model_class_key_map(VariantAlternate,
            variant_key_map[MAP_KEY__ALTERNATE])

    _update_model_class_key_map(VariantEvidence,
            variant_key_map[MAP_KEY__EVIDENCE])


@login_required
@require_POST
def modify_variant_in_set_membership(request):
    """Action that handles modifying the membership of a Variant in a
    VariantSet.
    """
    request_data = json.loads(request.body)

    # Make sure the required keys are present.
    # Validation.
    REQUIRED_KEYS = [
            'refGenomeUid',
            'variantSetAction',
            'variantSetUid']
    if not all(key in request_data for key in REQUIRED_KEYS):
        return HttpResponseBadRequest("Invalid request. Missing keys.")

    add_all_matching_filter = False
    if ('isAllMatchingFilterSelected' in request_data and
            request_data['isAllMatchingFilterSelected']):
        add_all_matching_filter = True
    else:
        if not 'variantUidList' in request_data:
            return HttpResponseBadRequest("Invalid request. Missing keys.")

    # Get the project and verify that the requesting user has the
    # right permissions.
    reference_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=request_data.get('refGenomeUid'))

    # Perform the update.
    if add_all_matching_filter:
        update_memberships_result = (
                update_variant_in_set_memberships__all_matching_filter(
                        reference_genome,
                        request_data.get('variantSetAction'),
                        request_data.get('variantSetUid'),
                        request_data.get('filterString'),
                        request_data.get('isMelted')))
    else:
        update_memberships_result = update_variant_in_set_memberships(
                reference_genome,
                request_data.get('variantUidList'),
                request_data.get('variantSetAction'),
                request_data.get('variantSetUid'))

    return HttpResponse(json.dumps(update_memberships_result))


@login_required
@require_GET
def get_variant_set_list(request):

    if 'refGenomeUid' in request.GET:
        # Parse the GET params.
        ref_genome_uid = request.GET.get('refGenomeUid')

        reference_genome = get_object_or_404(ReferenceGenome,
                project__owner=request.user.get_profile(),
                uid=ref_genome_uid)

        # Grab the VariantSet data.
        variant_set_list = VariantSet.objects.filter(
                reference_genome=reference_genome)

        response_data = {
            'variant_set_list_json': adapt_model_to_frontend(VariantSet,
                    obj_list=variant_set_list)
        }
        return HttpResponse(json.dumps(response_data),
                content_type='application/json')

    elif 'projectUid' in request.GET:
        project_uid = request.GET.get('projectUid')

        # Lookup the model and verify the owner is the user
        project = get_object_or_404(Project,
                owner=request.user.get_profile(),
                uid=project_uid)

        response_data = adapt_model_to_frontend(VariantSet,
                {'reference_genome__project':project})

        return HttpResponse(response_data,
                content_type='application/json')


@login_required
@require_GET
def get_samples(request):

    project_uid = request.GET.get('projectUid')

    # Lookup the model and verify the owner is the user
    project = get_object_or_404(Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    response_data = adapt_experiment_samples_to_frontend({'project': project})

    return HttpResponse(response_data,
            content_type='application/json')


@login_required
def get_gene_list(request):
    """Returns the Gene view data, showing Genes and aggregated counts.
    """
    ag_uid = request.GET.get('alignmentGroupUid')

    alignment_group = get_object_or_404(AlignmentGroup,
            reference_genome__project__owner=request.user.get_profile(),
            uid=ag_uid)

    gene_view_list = lookup_genes(alignment_group)

    response_data = {
        'geneList': adapt_gene_list_to_frontend(gene_view_list, alignment_group)
    }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


@login_required
@require_GET
def refresh_materialized_variant_table(request):
    """Updates the materialized variant table corresponding to the
    ReferenceGenome whose uid is provided in the GET params.
    """
    # DEBUG: Profiling refresh time.
    # import time
    # profiling_time_start = time.time()

    ref_genome_uid = request.GET.get('refGenomeUid')
    reference_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    # NOTE: Call create() for now. It may be possible to make this quicker
    # by calling refresh().
    mvmvm = MeltedVariantMaterializedViewManager(reference_genome)
    mvmvm.create()

    # print 'REFRESH TOOK', time.time() - profiling_time_start

    return HttpResponse('ok')


@require_GET
@login_required
def export_variants_as_csv(request):
    """Handles a request to download variants in .csv format.
    """
    ag_uid = request.GET.get('alignment_group_uid')
    alignment_group = get_object_or_404(AlignmentGroup,
            reference_genome__project__owner=request.user.get_profile(),
            uid=ag_uid)

    filter_string = request.GET.get('filter_string', '')

    response = StreamingHttpResponse(
            export_melted_variant_view(alignment_group, filter_string),
            content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="variants.csv"'
    return response


@login_required
@require_GET
def get_alignment_groups(request):
    """Get list of AlignmentGroups for the provided ReferenceGenome uid.

    If the request has a refGenomeUid, only return alignments for that
    individual reference genome.

    If the request has a projectUid, return all alignments for that
    project.

    TODO(gleb): Clarify in comments why we have these two cases.
    """
    if 'refGenomeUid' in request.GET:

        # Parse the GET params.
        ref_genome_uid = request.GET.get('refGenomeUid')

        # Lookup the model and verify the owner is hte user
        reference_genome = get_object_or_404(ReferenceGenome,
                project__owner=request.user.get_profile(),
                uid=ref_genome_uid)

        alignment_group_list = AlignmentGroup.objects.filter(
                reference_genome=reference_genome).order_by('label')

        response_data = [{
            'label': ag.label,
            'uid': ag.uid
        } for ag in alignment_group_list]

        return HttpResponse(json.dumps(response_data),
                content_type='application/json')

    elif 'projectUid' in request.GET:

        # Parse the GET params.
        project_uid = request.GET.get('projectUid')

        # Lookup the model and verify the owner is hte user
        project = get_object_or_404(Project,
                owner=request.user.get_profile(),
                uid=project_uid)

        alignment_group_list = AlignmentGroup.objects.filter(
                reference_genome__project=project).order_by('-start_time')

        response_data = adapt_model_to_frontend(AlignmentGroup,
                obj_list=alignment_group_list)

        # Add bit to indicate whether any AlignmentGroups are running.
        # NOTE: We do this wonky json.loads(), modify, json.dumps() becaause
        # adapt_model_to_frontend() has the suboptimal interface of returning
        # a json packaged object. It would be better to change this, but would
        # require making this change safely everywhere else, but since we are
        # lacking test coverage I'm not going to do that right now.
        response_data_dict = json.loads(response_data)
        response_data_dict['clientShouldRefresh'] = _are_any_alignments_running(
                alignment_group_list)
        response_data = json.dumps(response_data_dict)

        return HttpResponse(response_data,
                content_type='application/json')


def _are_any_alignments_running(alignment_group_list):
    """Determines whether any alignments in the list are running.
    """
    for ag in alignment_group_list:
        if ag.status in AlignmentGroup.PIPELINE_IS_RUNNING_STATUSES:
            return True
    return False


@login_required
@require_POST
def alignment_groups_delete(request):
    """Deletes AlignmentGroups.
    """
    request_data = json.loads(request.body)
    uid_list = request_data.get('uidList', [])
    if len(uid_list) == 0:
        raise Http404

    # First make sure all the samples belong to this user.
    to_delete = AlignmentGroup.objects.filter(
            reference_genome__project__owner=request.user.get_profile(),
            uid__in=uid_list)
    if not len(to_delete) == len(uid_list):
        raise Http404

    # Validation successful, delete.
    to_delete.delete()

    # Return success response.
    return HttpResponse(json.dumps({}), content_type='application/json')


@login_required
@require_GET
def is_materialized_view_valid(request):
    """Checks whether the materialized view is valid for this ReferenceGenome.
    """
    ref_genome_uid = request.GET.get('refGenomeUid')
    reference_genome = get_object_or_404(ReferenceGenome,
                project__owner=request.user.get_profile(),
                uid=ref_genome_uid)
    response_data = json.dumps({
        'isValid': reference_genome.is_materialized_variant_view_valid
    })
    return HttpResponse(response_data, content_type='application/json')


@login_required
@require_GET
def get_ref_genomes(request):
    """Get list of AlignmentGroups for the provided ReferenceGenome uid.
    """
    # Parse the GET params.
    project_uid = request.GET.get('projectUid')

    # Lookup the model and verify the owner is the user
    project = get_object_or_404(
            Project,
            owner=request.user.get_profile(),
            uid=project_uid)

    filters = {'project': project}

    response_data = adapt_model_to_frontend(ReferenceGenome, filters)

    return HttpResponse(response_data, content_type='application/json')


@login_required
@require_GET
def get_contigs(request):
    """Get list of Contigs for the provided Project uid.
    """
    # Parse the GET params.
    ref_genome_uid = request.GET.get('refGenomeUid')
    alignment_group_uid = request.GET.get('alignmentGroupUid')

    sample_to_align_query = ExperimentSampleToAlignment.objects.filter(
            alignment_group__uid=alignment_group_uid)

    filters = {
            'parent_reference_genome': ReferenceGenome.objects.get(
                    uid=ref_genome_uid),
            'experiment_sample_to_alignment__in': sample_to_align_query
    }

    response_data = adapt_model_to_frontend(Contig, filters)

    return HttpResponse(response_data, content_type='application/json')


@login_required
@require_GET
def contigs_has_insertion_location(request):
    """Returns {'has_insertion_location': True} if any of the contigs
    in the passed list of contig uids have insertion location data
    """
    # Parse the GET params.
    request_data = json.loads(request.GET['data'])
    contig_uid_list = request_data.get('contigUidList', None)

    # Search metadata of filtered contigs for non-null
    # insertion location keys
    has_insertion_location = False
    for c in Contig.objects.filter(uid__in=contig_uid_list):
        if c.metadata.get('insertion_sequence_endpoints', None):
            has_insertion_location = True
            break

    result = {'has_insertion_location': has_insertion_location}
    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def contigs_find_insertion_location(request):
    """Attempts to find the placement parameters for the contigs
    in the passedd contig uid list.  If unable to place a contig,
    the specific error message is included in the format (label, error_string)
    in a list under the key 'error' in the response dict
    """
    request_data = json.loads(request.body)
    contig_uid_list = request_data.get('contigUidList', [])

    result = {}
    for contig in Contig.objects.filter(uid__in=contig_uid_list):
        ref_genome = contig.parent_reference_genome
        contig_fasta = get_dataset_with_type(contig,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
        with open(contig_fasta) as fh:
            seqrecord = SeqIO.parse(fh, 'fasta').next()
            insertion_data = find_contig_insertion_site(
                    ref_genome, seqrecord)

            if 'error_string' in insertion_data:
                if result.get('error', False):
                    result['error'].append(
                            (contig.label, insertion_data['error_string']))
                else:
                    result['error'] = [
                            (contig.label, insertion_data['error_string'])]
            else:
                contig.metadata['insertion_sequence_endpoints'] = (
                        insertion_data['contig_cassette_start_pos'],
                        insertion_data['contig_cassette_end_pos'])
                contig.metadata['ref_insertion_pos'] = insertion_data[
                        'ref_insertion_pos']
                contig.metadata['chromosome'] = insertion_data[
                        'ref_chromosome_seqrecord_id']
                contig.save()

    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_POST
def contigs_place_in_ref(request):
    """Incorporates the passed contigs into the reference they belong to
    to make a new version of the reference with the passed label. For now
    only incorporation of single contigs is supported
    """
    request_data = json.loads(request.body)
    contig_uid_list = request_data.get('contigUidList', [])
    new_genome_label = request_data.get('newGenomeLabel', '')

    result = {}
    for contig in Contig.objects.filter(uid__in=contig_uid_list):
        ref_genome = contig.parent_reference_genome

        start, end = contig.metadata.get('insertion_sequence_endpoints')
        placement_position_params = {
            'ref_insertion_pos': contig.metadata.get('ref_insertion_pos'),
            'ref_chromosome_seqrecord_id': contig.metadata.get('chromosome'),
            'contig_cassette_start_pos': start,
            'contig_cassette_end_pos': end
        }

        new_reference_genome_params = {
            'label': new_genome_label
        }

        contig_fasta = get_dataset_with_type(contig,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
        with open(contig_fasta) as fh:
            seqrecord = SeqIO.parse(fh, 'fasta').next()
            place_cassette(ref_genome, seqrecord,
                    placement_position_params, new_reference_genome_params)

    return HttpResponse(json.dumps(result), content_type='application/json')


@login_required
@require_GET
def get_single_ref_genome(request):
    reference_genome_uid = request.GET.get('referenceGenomeUid')
    response_data = adapt_model_to_frontend(
        Chromosome, {'reference_genome__uid': reference_genome_uid})

    return HttpResponse(response_data, content_type='application/json')


@login_required
@require_POST
def create_variant_set(request):
    # Get the params.
    ref_genome_uid = request.POST.get('refGenomeUid', '')
    variant_set_name = request.POST.get('variantSetName', '')
    create_set_type = request.POST.get('createSetType', '')

    # Basic validation.
    try:
        assert create_set_type in ['from-file', 'empty']
        assert ref_genome_uid != '', "Must provide Reference Genome"
        assert variant_set_name != '', "Must provide Variant Set name"
    except AssertionError as e:
        return HttpResponseBadRequest(str(e))

    # Model lookup / validation.
    ref_genome = get_object_or_404(
            ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    # Create new variant set, depending on type of form submitted.
    if create_set_type == 'from-file':
        result = _create_variant_set_from_file(
                request, ref_genome, variant_set_name)
    else:
        result = _create_variant_set_empty(ref_genome, variant_set_name)

    return HttpResponse(json.dumps(result), content_type='application/json')


def _create_variant_set_from_file(request, ref_genome, variant_set_name):
    """Creates a variant set from uploaded vcf file.

    Returns:
        Dictionary with keys:
            * error_str: Either empty string or description of error that occurred
    """
    error_string = ''

    path = default_storage.save('tmp/tmp_varset.vcf',
            ContentFile(request.FILES['vcfFile'].read()))
    variant_set_file = os.path.join(settings.MEDIA_ROOT, path)

    try:
        file_variant_set = import_variant_set_from_vcf(ref_genome, variant_set_name,
                variant_set_file)
    except Exception as e:
        error_string = 'Import error: ' + str(e)
    finally:
        os.remove(variant_set_file)

    result = {
        'error': error_string,
    }

    return result


def _create_variant_set_empty(ref_genome, variant_set_name):
    """Creates an empty variant set.

    A VariantSet with the given name can't exist already.

    Returns:
        Dictionary with keys:
            * error_str: Either empty string or description of error that occurred
            * variant_set_uid: uid of the new VariantSet
    """
    exists_set_with_same_name = bool(VariantSet.objects.filter(
        reference_genome=ref_genome,
        label=variant_set_name).count())

    if exists_set_with_same_name:
        error_string = 'Variant set %s exists' % variant_set_name
        result = {
            'error': error_string,
        }
    else:
        error_string = ''
        empty_variant_set = VariantSet.objects.create(
                reference_genome=ref_genome,
                label=variant_set_name)
        result = {
            'error': error_string,
            'variantSetUid': empty_variant_set.uid
        }

    return result


@login_required
def print_mage_oligos_for_variant_set(request):
    variant_set_uid = request.GET.get('variantSetUid')
    variant_set = get_object_or_404(VariantSet,
            reference_genome__project__owner=request.user.get_profile(),
            uid=variant_set_uid)
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="oligos.csv"'
    repliation_origin_params = ReplicationOriginParams(
            request.GET.get('repOriginStart'),
            request.GET.get('repOriginEnd'),
            request.GET.get('repTerminusStart'),
            request.GET.get('repTerminusEnd'))
    print_mage_oligos(variant_set, response, 'o_', repliation_origin_params,
            experiment_dir=request.GET.get('experimentDir'))
    return response


@login_required
@require_POST
def generate_new_ref_genome_for_variant_set(request):
    variant_set_uid = request.POST.get('variantSetUid')
    variant_set = get_object_or_404(VariantSet,
            reference_genome__project__owner=request.user.get_profile(),
            uid=variant_set_uid)
    new_ref_genome_label = request.POST.get('refGenomeLabel')
    ref_genome_maker_params = {
        'label': new_ref_genome_label
    }

    error_string = ''
    try:
        new_ref_genome = generate_new_reference_genome(
                variant_set, ref_genome_maker_params)
    except ValidationException as e:
        error_string = str(e)

    if not error_string:
        assert new_ref_genome
        result = {
            'redirect': reverse(
                    'main.views.reference_genome_view',
                    args=(new_ref_genome.project.uid, new_ref_genome.uid)),
        }
    else:
        result = {
            'error': error_string
        }
    return HttpResponse(json.dumps(result), content_type='application/json')


@require_GET
@login_required
def generate_contigs(request):
    """
    Generates and begins download of a fasta file of contigs assembled from
    unmapped and split reads of the passed ExperimentSampleToAlignment
    """

    # Retrieve ExperimentSampleToAlignment
    sample_alignment_uid = request.GET.get('sampleAlignmentUid')
    experiment_sample_to_alignment = get_object_or_404(
            ExperimentSampleToAlignment,
            alignment_group__reference_genome__project__owner=(
                    request.user.get_profile()),
            uid=sample_alignment_uid)

    # Get reference genome
    reference_genome = (
            experiment_sample_to_alignment.alignment_group.reference_genome)

    # Generate name for contigs
    sample_label = experiment_sample_to_alignment.experiment_sample.label
    ref_label = reference_genome.label
    contig_label_base = '_'.join(
            [ref_label, sample_label])

    # Generate a list of fasta file paths to the contigs
    contig_filepaths = assembly.generate_contigs(
            experiment_sample_to_alignment, contig_label_base)

    # Check if contigs exist
    are_no_contigs = all([os.stat(contig_filepath).st_size == 0
            for contig_filepath in contig_filepaths])

    result = {'is_contig_file_empty': are_no_contigs}

    return HttpResponse(
        json.dumps(result), content_type='application/json')


if settings.S3_ENABLED:
    @login_required
    def import_reference_genome_s3(request, project_uid):
        if request.method == 'POST':
            project = get_object_or_404(Project, owner=request.user.get_profile(),
                uid=project_uid)
            s3file_id = request.POST['s3file_id']
            s3file = S3File.objects.get(pk=s3file_id)
            import_reference_genome_from_s3(
                        project,
                        request.POST['refGenomeLabel'],
                        s3file,
                        request.POST['importFileFormat'])
            
        return HttpResponse("", content_type='text/plain')

    @login_required
    @require_POST
    def parse_targets_file_s3(request, project_uid):
        project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
        s3file_id = request.POST['s3file_id']
        s3file = S3File.objects.get(pk=s3file_id)

        csv_data = s3_get_string(s3file.key)
        csv_io = StringIO(csv_data)
        sample_filenames = []

        try:
            valid_rows = parse_targets_file(csv_io, remove_directory_path=True)
            for field_name, field_value in valid_rows.iteritems():
                if 'Path' in field_name:
                    sample_filenames.append(field_value)
        except AssertionError as e:
            return HttpResponse(json.dumps({
                    'error': str(e)
                }), content_type='application/json')
        except:
            import traceback
            return HttpResponse(json.dumps({
                    'error': traceback.format_exc()
                }), content_type='application/json')

        if len(list(set(sample_filenames))) != len(sample_filenames):
            return HttpResponse(json.dumps({
                    'error': "Targets file contains sample files with same names."
                }), content_type='application/json')

        return HttpResponse(json.dumps({
                'targets_file_rows': valid_rows,
                'sample_filenames': sample_filenames
            }), content_type='application/json')

    @login_required
    @require_POST
    def process_sample_files_s3(request, project_uid):
        project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
        data = json.loads(request.raw_post_data)
        s3files = []
        for f in data['sample_files'].values():
            s3files.append(S3File.objects.get(pk=int(f['sid'])))

        import_samples_from_s3(project, data['targets_file_rows'], s3files)

        return HttpResponse(json.dumps({
                'targets_file_rows': data['targets_file_rows'],
                'sample_files': data['sample_files']
            }), content_type='application/json')
