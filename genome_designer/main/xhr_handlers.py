"""
Methods that handle Ajax requests from the frontend.

This module was created in response to views.py getting quite big, and so a
reasonable separation point is to separate page actions from Ajax actions.
"""

import copy
import csv
import json
from StringIO import StringIO
import time

from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.http import Http404
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_http_methods
from django.views.decorators.http import require_GET, require_POST

from main.adapters import adapt_model_or_modelview_list_to_frontend
from main.adapters import adapt_model_to_frontend
from main.data_util import lookup_variants
from main.model_views import adapt_variant_to_frontend
from main.model_views import GeneView
from main.models import AlignmentGroup
from main.models import Project
from main.models import ReferenceGenome
from main.models import Region
from main.models import VariantCallerCommonData
from main.models import VariantAlternate
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import S3File
from scripts.data_export_util import export_melted_variant_view
from scripts.dynamic_snp_filter_key_map import MAP_KEY__COMMON_DATA
from scripts.dynamic_snp_filter_key_map import MAP_KEY__ALTERNATE
from scripts.dynamic_snp_filter_key_map import MAP_KEY__EVIDENCE
from variants.common import extract_filter_keys
from variants.materialized_variant_filter import get_variants_that_pass_filter
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager
from variants.variant_sets import update_variant_in_set_memberships

if settings.S3_ENABLED:
    from scripts.import_util import parse_targets_file, import_reference_genome_from_s3, import_samples_from_s3
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
            'alt': variant.get_variants_as_string(),
        })
    return response


# Key in the GET params containing the string for filtering the variants.
VARIANT_LIST_REQUEST_KEY__FILTER_STRING = 'variantFilterString'
VARIANT_LIST_REQUEST_KEY__PROJECT_UID = 'projectUid'
VARIANT_LIST_REQUEST_KEY__REF_GENOME_UID = 'refGenomeUid'

VARIANT_LIST_RESPONSE_KEY__LIST = 'variant_list_json'
VARIANT_LIST_RESPONSE_KEY__TOTAL = 'num_total_variants'
VARIANT_LIST_RESPONSE_KEY__SET_LIST = 'variant_set_list_json'
VARIANT_LIST_RESPONSE_KEY__KEY_MAP = 'variant_key_filter_map_json'
VARIANT_LIST_RESPONSE_KEY__ERROR = 'error'


@login_required
@require_GET
def get_variant_list(request):
    """Returns a list of Variants, filtered by any filter parameters contained
    in the request.
    """
    # Parse the GET params.
    ref_genome_uid = request.GET.get('refGenomeUid')
    project_uid = request.GET.get('projectUid')

    # Get model and verify permisssions.
    reference_genome = get_object_or_404(ReferenceGenome,
            project__uid=project_uid, uid=ref_genome_uid)

    # Dictionary to hold all query specific parameters.
    query_args = {}

    # Get inputs to perform the query for Variants data.
    # TODO: Combine with saved filter string.
    query_args['filter_string'] = request.GET.get(
            VARIANT_LIST_REQUEST_KEY__FILTER_STRING, '')

    # Determine whether melted or cast view.
    query_args['is_melted'] = request.GET.get('melt', 0) == '1'

    # Get optional column to sort by.
    query_args['sort_by_column'] = request.GET.get('sortBy', '')

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
        # Determine the visible keys.
        visible_key_names = _determine_visible_field_names(request,
                query_args['filter_string'], reference_genome)

        # Get the list of Variants (or melted representation) to display.
        lookup_variant_result = lookup_variants(query_args, reference_genome)
        variant_list = lookup_variant_result.result_list
        num_total_variants = lookup_variant_result.num_total_variants

        # Adapt the Variants to display for the frontend.
        variant_list_json = adapt_variant_to_frontend(variant_list,
                reference_genome, visible_key_names,
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
                variant_key_map_with_active_fields_marked, visible_key_names)

        # Package up the response.
        response_data = {
            VARIANT_LIST_RESPONSE_KEY__LIST: variant_list_json,
            VARIANT_LIST_RESPONSE_KEY__TOTAL: num_total_variants,
            VARIANT_LIST_RESPONSE_KEY__SET_LIST: adapt_model_to_frontend(VariantSet,
                    obj_list=variant_set_list),
            VARIANT_LIST_RESPONSE_KEY__KEY_MAP: json.dumps(
                    variant_key_map_with_active_fields_marked)
        }

    except Exception as e:
        # TODO: More readable error reporting.
        exception_as_string = str(type(e)) + ' ' + str(e)
        response_data = {
            VARIANT_LIST_RESPONSE_KEY__ERROR: exception_as_string
        }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS = 'visibleKeyNames'


def _determine_visible_field_names(request, filter_string, ref_genome):
    """Determine which fields to show.
    """
    # Get visible keys explicitly marked in the UI by the user.
    if VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS in request.GET:
        visible_key_names = json.loads(request.GET.get(
                VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS))
    else:
        visible_key_names = []

    # Also show keys in the filter string.
    fields_from_filter_string = extract_filter_keys(filter_string, ref_genome)

    return list(set(visible_key_names) | set(fields_from_filter_string))


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
@require_http_methods(['POST'])
def modify_variant_in_set_membership(request):
    """Action that handles modifying the membership of a Variant in a
    VariantSet.
    """
    request_data = json.loads(request.body)

    # Make sure the required keys are present.
    REQUIRED_KEYS = [
            'refGenomeUid',
            'projectUid',
            'variantUidList',
            'variantSetAction',
            'variantSetUid']

    # Validate the request.
    if not all(key in request_data for key in REQUIRED_KEYS):
        return HttpResponseBadRequest("Invalid request. Missing keys.")

    ref_genome_uid = request_data.get('refGenomeUid')
    project_uid = request_data.get('projectUid')

    # Get the project and verify that the requesting user has the
    # right permissions.
    project = get_object_or_404(Project, owner=request.user.get_profile(),
            uid=project_uid)
    reference_genome = get_object_or_404(ReferenceGenome, project=project,
            uid=ref_genome_uid)

    # Add or remove the variants to the set, as per variantSetAction.
    update_memberships_result = update_variant_in_set_memberships(
            reference_genome,
            request_data.get('variantUidList'),
            request_data.get('variantSetAction'),
            request_data.get('variantSetUid'))

    return HttpResponse(json.dumps(update_memberships_result))


@login_required
@require_GET
def get_variant_set_list(request):
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
                obj_list=variant_set_list),
    }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


@login_required
def get_gene_list(request):
    ref_genome_uid = request.GET.get('refGenomeUid')

    reference_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    region_list = Region.objects.filter(
            reference_genome=reference_genome,
            type=Region.TYPE.GENE)
    gene_view_list = [GeneView(region) for region in region_list]

    response_data = {
        'geneList': adapt_model_or_modelview_list_to_frontend(gene_view_list)
    }

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


@login_required
@require_GET
def refresh_materialized_variant_table(request):
    """Updates the materialized variant table corresponding to the
    ReferenceGenome whose uid is provided in the GET params.
    """
    profiling_time_start = time.time()

    ref_genome_uid = request.GET.get('refGenomeUid')
    reference_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    # NOTE: Call create() for now. It may be possible to make this quicker
    # by calling refresh().
    mvmvm = MeltedVariantMaterializedViewManager(reference_genome)
    mvmvm.create()

    print 'REFRESH TOOK', time.time() - profiling_time_start

    return HttpResponse('ok')


@login_required
def export_variants_as_csv(request):
    """Handles a request to download variants in .csv format.
    """
    ref_genome_uid = request.GET.get('ref_genome_uid')
    reference_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    # NOTE: Currently a no-op.
    variant_id_list = []

    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="variants.csv"'
    export_melted_variant_view(reference_genome, variant_id_list, response)
    return response


@login_required
@require_GET
def get_alignment_groups_for_ref_genome(request):
    """Get list of AlignmentGroups for the provided ReferenceGenome uid.
    """
    # Parse the GET params.
    ref_genome_uid = request.GET.get('refGenomeUid')

    # Lookup the model and verify the owner is hte user
    reference_genome = get_object_or_404(ReferenceGenome,
            project__owner=request.user.get_profile(),
            uid=ref_genome_uid)

    alignment_group_list = AlignmentGroup.objects.filter(
            reference_genome=reference_genome)
    response_data = [{
        'label': ag.label,
        'uid': ag.uid
    } for ag in alignment_group_list]

    return HttpResponse(json.dumps(response_data),
            content_type='application/json')


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
