"""
Classes that describe how a particular model should be viewed.
"""

from collections import defaultdict, OrderedDict
import json
import string

from django.core.urlresolvers import reverse

from main.models import AlignmentGroup
from main.models import Variant
from main.model_view_utils import create_variant_links_field
from main.model_view_utils import get_jbrowse_track_names
from main.model_view_utils import create_alt_flag_field
from utils import titlecase_spaces
from variants.common import generate_key_to_materialized_view_parent_col
from variants.common import validate_key_against_map
from variants.melted_variant_schema import CAST_SCHEMA_KEY__TOTAL_SAMPLE_COUNT
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__ALT
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__HET
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__CHROMOSOME
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__ES_LABEL
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__ES_UID
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__POSITION
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__REF
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__UID
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VA_ID
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VS_LABEL
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__VS_UID


ALL_VS_LABEL_KEY = 'all_variant_set_label'
ALL_VS_UID_KEY = 'all_variant_set_uid'


def adapt_non_recursive(obj_list, field_dict_list, reference_genome, melted):
    """Adapts of list of objects that doesn't require recursive calling.

    Returns:
        JSON string representation of frontend objects.
    """
    # We want a list of all VCF tracks for jbrowse. Makes sense to do it
    # once and then pass the strings to each variant object for the frontend.
    # This is a little hacky, but whatever.
    jbrowse_track_names = get_jbrowse_track_names(reference_genome)

    # HACK: Only one AlignmentGroup right now.
    associated_alignment_groups = AlignmentGroup.objects.filter(
            reference_genome=reference_genome)
    if len(associated_alignment_groups) > 0:
        hack_single_alignment_group = associated_alignment_groups[0]
    else:
        hack_single_alignment_group = None

    # Aggregate list of objects that are ready for display by the frontend.
    fe_obj_list = []
    for melted_variant_obj in obj_list:
        visible_field_pairs = []

        # If there is an empty row (no ExperimentSample associated),
        # then all the counts will be off by one, so we need to decrement
        # them.
        if not melted:
            if None in melted_variant_obj[MELTED_SCHEMA_KEY__ES_UID]:
                maybe_dec = 1
            else:
                maybe_dec = 0
            assert maybe_dec >= 0, "maybe_dec should be positive"
        else:
            maybe_dec = 0

        for fdict in field_dict_list:
            value = None
            field = fdict['field']

            # HACK: Special handling for certain fields.
            if field == 'links':
                value = create_variant_links_field(
                        melted_variant_obj, reference_genome,
                        jbrowse_track_names, hack_single_alignment_group)
            elif field == MELTED_SCHEMA_KEY__VS_LABEL:
                value = _adapt_variant_set_label_field(
                        melted_variant_obj, reference_genome.project.uid,
                        melted, reference_genome, hack_single_alignment_group)
            elif field == MELTED_SCHEMA_KEY__REF and not melted:
                value = (melted_variant_obj[MELTED_SCHEMA_KEY__REF] + ' (%d)' %
                        melted_variant_obj[MELTED_SCHEMA_KEY__ALT].count(None))
            elif field == MELTED_SCHEMA_KEY__ALT:
                value = create_alt_flag_field(melted_variant_obj, melted, maybe_dec)
            elif field == CAST_SCHEMA_KEY__TOTAL_SAMPLE_COUNT:
                value = (melted_variant_obj[
                        CAST_SCHEMA_KEY__TOTAL_SAMPLE_COUNT] - maybe_dec)
            elif fdict.get('is_subkey', False):
                assert 'parent_col' in fdict
                parent_dict_or_list = melted_variant_obj.get(fdict['parent_col'].upper(), {})
                if isinstance(parent_dict_or_list, dict):
                    # melted
                    parent_dict = parent_dict_or_list
                    value = adapt_melted_object_field(
                            parent_dict.get(field, ''), fdict)

                elif isinstance(parent_dict_or_list, list):
                    # cast
                    value = adapt_cast_object_list_field(
                            parent_dict_or_list, fdict)
            else:
                value = melted_variant_obj.get(field, '')
            # Pass empty string when no value present.
            if value is None:
                value = ''

            # Append (key, value) to list so we can create dict at the end.
            # Order doesn't matter here since it's determined in the
            # field config object constructed separately below.
            visible_field_pairs.append((field, value))
        # Append catch all INFO column field
        # TODO use variables instead of magic strings
        va_data = melted_variant_obj['VA_DATA']

        # MAJOR HACK ALERT: We need to rethink the whole thing.
        if isinstance(va_data, list):
            # Default.
            representative_va = va_data[0]

            # Try to do better. First one might not have data.
            for va in va_data:
                if not va:
                    continue
                if va.get('INFO_EFF_AA', '') or va.get('INFO_SVTYPE', ''):
                    representative_va = va
                    break

            # Use this as va_data from here on.
            va_data = representative_va

        if va_data:
            if 'INFO_SVTYPE' in va_data:
                # is SV: make info of the form "SV [type] [length]"
                key = 'INFO'
                value = 'SV {svtype} {svlen}'.format(
                        svtype=va_data['INFO_SVTYPE'],
                        svlen=va_data.get('INFO_SVLEN', ''))
                visible_field_pairs.append((key, value))
            else:
                # is SNP: make info equal to the AA field
                visible_field_pairs.append(('INFO', va_data.get('INFO_EFF_AA', '')))
        else:
            # unknown: just leave info field blank
            visible_field_pairs.append(('INFO', ''))

        fe_obj_list.append(dict(visible_field_pairs))

    # Create the config dict that tells DataTables js how to display each col.
    obj_field_config = []
    # save idxes of fields we want to be last
    last_idxes = []
    # Keep track of field indexes for moving the 'last' fields
    # need this counter instead of enumerate() to prevent
    # counting hidden fields
    displayed_idx = 0
    for fdict in field_dict_list:
        if fdict.get('hide', False):
            continue
        config_obj = {
            'mData': fdict['field'],
            'sTitle': fdict.get('verbose',
                    string.capwords(fdict['field'], '_').replace('_', ' ')),
            'bSortable': fdict.get('sortable', False)
        }
        if fdict.get('sWidth'):
            config_obj['sWidth'] = fdict.get('sWidth')
        obj_field_config.append(config_obj)
        if fdict.get('last'):
            last_idxes.append(displayed_idx)
        displayed_idx += 1

    # Add metadata for INFO column
    obj_field_config.append({
      'mData': 'INFO',
      'sTitle': 'Info',
      'bSortable': False
      })
    displayed_idx += 1

    # Finally, move these fields to the end.
    for i in last_idxes:
        field_i = obj_field_config.pop(i)
        obj_field_config.append(field_i)

    return json.dumps({
        'obj_list': fe_obj_list,
        'field_config': obj_field_config
    })

def _create_label_for_variant_object(variant_as_dict):
    # Generate label from variant data.
    # Using ordered dicts with None as an ordered set
    # NOTE: Not using this currently, but we might want it in the future.

    position = str(variant_as_dict[MELTED_SCHEMA_KEY__POSITION])

    try:
        va_data = variant_as_dict['VA_DATA']

        genes = OrderedDict()
        effs = OrderedDict()
        aas = OrderedDict()
        impacts = OrderedDict()
        has_effect = False

        # Cast view
        if isinstance(va_data, list):
            for alt in va_data:

                if alt is None: continue

                if 'INFO_EFF_EFFECT' in alt:
                    has_effect = True

                gene = alt.get('INFO_EFF_GENE', None)

                if not gene: continue

                eff = string.lower(alt['INFO_EFF_EFFECT'].replace('_',' '))
                impact = string.lower(alt['INFO_EFF_IMPACT'].replace('_',' '))
                aa = alt['INFO_EFF_AA']
                if gene is not 'None': genes[gene] = None
                if aa is not 'None': aas[aa] = None
                if eff is not 'None': effs[eff] = None
                if impact is not 'None': impacts[impact] = None

        # Melted view
        else:
            if va_data is not None:

                if 'INFO_EFF_EFFECT' in va_data:
                        has_effect = True

                if va_data.get('INFO_EFF_GENE', None) is not None:

                    genes[va_data['INFO_EFF_GENE']] = None
                    aas[va_data['INFO_EFF_AA']] = None
                    effs[(string.lower(
                            va_data['INFO_EFF_EFFECT'].replace('_',' ')))] = None
                    impacts[(string.lower(
                            va_data['INFO_EFF_IMPACT'].replace('_',' ')))] = None

        # Generate the label.
        if len(genes) == 0:
            if va_data and has_effect:
                label = label = position + ': intergenic'
            else:
                label = label = position + ': unannotated/ref'
        else:
            label =  '{:s}:{:s}, {:s} ({:s})'.format(
                    '/'.join(genes.keys()),
                    '/'.join(aas.keys()),
                    ','.join(effs.keys()),
                    ','.join(impacts.keys()))

    except Exception as e:
        print str(e)
        print va_data
        label = position + ': error'

    return label


def _adapt_variant_set_label_field(variant_as_dict, project_uid, melted,
        reference_genome, alignment_group):
    """Constructs the labels as anchors that link to the single variant view.
    """
    if melted:
        return _adapt_variant_set_label_field__melted(
                variant_as_dict, project_uid, reference_genome,
                alignment_group)
    else:
        return _adapt_variant_set_label_field__cast(
                variant_as_dict, project_uid, reference_genome,
                alignment_group)


def _adapt_variant_set_label_field__melted(variant_as_dict, project_uid,
        reference_genome, alignment_group=None):
    # Build a dictionary of individual HTML string anchors mapped by label,
    # so we can sort it at the very end.
    variant_set_anchor_map = {}

    # Also get a map of uid to field.
    uid_to_label_map = dict(zip(
            variant_as_dict[ALL_VS_UID_KEY],
            variant_as_dict[ALL_VS_LABEL_KEY]))
    for uid, label in uid_to_label_map.items():
        if uid is None or label is None:
            continue

        # This is the link to the variant set view.
        variant_set_href = _create_variant_set_analyze_view_link(
                project_uid, uid, reference_genome, alignment_group)

        # If the variant set is for this sample, then it will be filled,
        # Otherwise, it will be outlined. Cast view always uses outline.
        variant_set_classes = ['gd-variant-set-badge']
        if not label in variant_as_dict[MELTED_SCHEMA_KEY__VS_LABEL]:
            variant_set_classes.append('outline')

        # Build the html to display.
        variant_set_html = (
                "<a href='" + variant_set_href + "'" +
                "class='" + ' '.join(variant_set_classes) + "'>" +
                label + "</a>")

        # Store it in the map.
        variant_set_anchor_map[label] = variant_set_html

    # Finally, return these anchors alphabetically sorted by label.
    sorted_set_labels = sorted(variant_set_anchor_map.keys())
    return ' '.join(
            [variant_set_anchor_map[i] for i in sorted_set_labels])


def _adapt_variant_set_label_field__cast(variant_as_dict, project_uid,
        reference_genome, alignment_group=None):
    # If there is an empty row (no ExperimentSample associated),
    # then all the counts will be off by one, so we need to decrement
    # them.
    if None in variant_as_dict[MELTED_SCHEMA_KEY__ES_UID]:
        maybe_dec = 1
    else:
        maybe_dec = 0
    assert maybe_dec >= 0, "Remaining code expects maybe_dec to be positive"

    # Bucket the labels and count occurrences.
    variant_set_label_to_count_map = defaultdict(lambda: 0)
    for label in variant_as_dict[MELTED_SCHEMA_KEY__VS_LABEL]:
        if not label:
            continue
        variant_set_label_to_count_map[label] += 1

    # Apply decrement if applicable.
    for label, count in variant_set_label_to_count_map.iteritems():
        variant_set_label_to_count_map[label] = count - maybe_dec

    # Build a map from label to uid so we can make hrefs.
    variant_set_label_to_uid_map = dict(zip(
            variant_as_dict[MELTED_SCHEMA_KEY__VS_LABEL],
            variant_as_dict[MELTED_SCHEMA_KEY__VS_UID]))

    # Build a dictionary of individual HTML string anchors mapped by label.
    # Sort at the end.
    variant_set_anchor_map = {}
    for label, count in variant_set_label_to_count_map.iteritems():
        uid = variant_set_label_to_uid_map[label]

        # This is the link to the variant set view.
        variant_set_href = _create_variant_set_analyze_view_link(
                project_uid, uid, reference_genome, alignment_group)

        # If the variant set is for this sample, then it will be filled,
        # Otherwise, it will be outlined. Cast view always uses outline.
        variant_set_classes = ['gd-variant-set-badge', 'outline', ' gd-nowrap']

        # If there are no samples associated, then don't display ': 0'
        if count > 0:
            text_content = label + ": " + str(count)
        else:
            text_content = label

        # Build the html to display.
        variant_set_html = (
                "<a href='" + variant_set_href + "'" +
                "class='" + ' '.join(variant_set_classes) + "'>" +
                text_content + "</a>")

        # Store it in the map.
        variant_set_anchor_map[label] = variant_set_html

    # Finally, return these anchors alphabetically sorted by label.
    sorted_set_labels = sorted(variant_set_anchor_map.keys())
    return ' '.join(
            [variant_set_anchor_map[i] for i in sorted_set_labels])


def _create_variant_set_analyze_view_link(project_uid, variant_set_uid,
        reference_genome, alignment_group=None):
    """Create link to Analyze view filtered by this variant set.
    """
    reverse_args = [reference_genome.project.uid]
    if alignment_group is not None:
        reverse_args += [alignment_group.uid, 'variants']
    root_href = reverse('main.views.tab_root_analyze', args=reverse_args)
    filter_part = '?filter=VARIANT_SET_UID=%s&melt=0' % (variant_set_uid,)
    return root_href + filter_part


# Field that provides links to melted view, JBrowse, etc.
LINKS_FIELD = {
    'field': 'links',
    'sWidth': '48px'
}

MELTED_VARIANT_FIELD_DICT_LIST = [
    LINKS_FIELD,
    {'field': MELTED_SCHEMA_KEY__ES_LABEL, 'verbose': 'Sample'},
    {'field': MELTED_SCHEMA_KEY__CHROMOSOME},
    {'field': MELTED_SCHEMA_KEY__POSITION},
    {'field': MELTED_SCHEMA_KEY__REF},
    {'field': MELTED_SCHEMA_KEY__ALT},
    {'field': MELTED_SCHEMA_KEY__HET, 'hide': True},
    {'field': MELTED_SCHEMA_KEY__VS_LABEL, 'verbose': 'Sets', 'last': True},
    {'field': MELTED_SCHEMA_KEY__ES_UID, 'hide': True},
    {'field': MELTED_SCHEMA_KEY__UID, 'hide': True}
]

CAST_VARIANT_FIELD_DICT_LIST = [
    LINKS_FIELD,
    {'field': MELTED_SCHEMA_KEY__CHROMOSOME},
    {
        'field': MELTED_SCHEMA_KEY__POSITION,
        'sortable': True
    },
    {'field': MELTED_SCHEMA_KEY__REF},
    {'field': MELTED_SCHEMA_KEY__ALT},
    {'field': MELTED_SCHEMA_KEY__HET, 'hide': True},
    {
        'field': CAST_SCHEMA_KEY__TOTAL_SAMPLE_COUNT,
        'verbose': '# Samples',
        'sortable': True
    },
    {'field': MELTED_SCHEMA_KEY__VS_LABEL, 'verbose': 'Sets', 'last': True},
    {'field': MELTED_SCHEMA_KEY__UID, 'hide': True}
]

# Fields that are added if they are available for the data.
# NOTE: The source table for these must be included in
# MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS.
OPTIONAL_DEFAULT_FIELDS = [
    {'field': 'INFO_EFF_GENE', 'verbose': 'Gene', 'format': 'gather'}, # va_data
    {'field': 'INFO_EFF_IMPACT', 'verbose': 'Impact', 'format': 'gather',
            'recase':'title'}, # va_data
]


def get_all_fields(reference_genome, visible_key_names, melted=False):
    """Gets the list of columns that will be displayed on DataTables, in
    display order.

    Args:
        reference_genome: The ReferenceGenome. We pass this so that we have a
            handle to the valid key map so we can handle additional
            visible key names properly.
        visible_key_names: List of key names that we want to include in the
            response.

    Returns:
        Tuple of default fields, optional fields, and additional fields,
        in display order.
    """
    if melted:
        default_field_dict_list = MELTED_VARIANT_FIELD_DICT_LIST
    else:
        default_field_dict_list = CAST_VARIANT_FIELD_DICT_LIST

    # Some of the keys are actually inside of catch-all data objects. Our
    # materialized view has separate columns for each of these, so we need
    # a structure that maps from key name to column name.
    # NOTE: We enforce unique keys when dynamically generating
    # ReferenceGenome.variant_key_map.
    key_to_parent_map = generate_key_to_materialized_view_parent_col(
            reference_genome)

    # Maybe add additional optional fields.
    optional_default_field_dict_list = []
    for field_dict in OPTIONAL_DEFAULT_FIELDS:
        field = field_dict['field']
        if validate_key_against_map(field_dict['field'],
                reference_genome.variant_key_map):
            new_field_dict = dict(field_dict.items() +
                    _prepare_visible_key_name_for_adapting_to_fe(
                            field, key_to_parent_map).items())
            optional_default_field_dict_list.append(new_field_dict)

    # Prepare additional fields for adaptation.
    additional_visible_field_dict_list = _prepare_additional_visible_keys(
            visible_key_names, key_to_parent_map)

    combined_list = (default_field_dict_list +
            optional_default_field_dict_list +
            additional_visible_field_dict_list)

    # Before we return, put last fields last.
    final_list = []
    last_fields = []
    for field in combined_list:
        if 'last' in field:
            last_fields.append(field)
        else:
            final_list.append(field)
    return final_list + last_fields


def adapt_variant_to_frontend(obj_list, reference_genome, visible_key_names,
        melted=False):
    """Convert the dictionary objects returned by Postgres into the form that
    the Datatables.js component can display.

    Args:
        obj_list: List of dictionaries representing melted entities returned
            by a raw SQL call directly to Postgres.
        reference_genome: The ReferenceGenome. We pass this so that we have a
            handle to the valid key map so we can handle additional
            visible key names properly.
        visible_key_names: List of key names that we want to include in the
            response.

    Returns:
        JSON string with the objects in the form that can be drawn by the
        Datatables component.
    """
    if melted:
        modified_obj_list = _modify_obj_list_for_variant_set_display(obj_list)
    else:
        modified_obj_list = obj_list

    all_field_dict_list = get_all_fields(
            reference_genome, visible_key_names, melted)
    
    return adapt_non_recursive(modified_obj_list, all_field_dict_list,
            reference_genome, melted)

def _modify_obj_list_for_variant_set_display(obj_list):
    # Before we adapt the fields, we need to do some special handling because
    # of the way that we structure the materialized view for the melted variant
    # data. Specifically, we sometimes include an extra row for each Variant,
    # which includes VariantSets that are not associated with any
    # ExperimentSample. (See the part of the materialized view build query that
    # follows UNION). The only case in which we want to show this catch-all
    # row is whe there are no samples associated with this variant at all.
    # Otherwise, we just want to grab the extra VariantSet data and show it
    # in the rows that do have samples associated, with a UI affordance that
    # indicates that the VariantSet is not explicitly associated with that
    # Sample, even though it is associated with that Variant.

    # The steps to implementing this strategy are:
    #     1. Triage how many rows there are for each Variant.
    #     2. For each redundant row, get rid of it, but save a reference to
    #        its variant_set_labels.

    
    def _make_catch_all_key(obj):
        return (obj[MELTED_SCHEMA_KEY__UID], obj[MELTED_SCHEMA_KEY__VA_ID])

    # First count rows for each Variant.
    variant_uid_to_count_dict = defaultdict(lambda: 0)
    for obj in obj_list:
        key = _make_catch_all_key(obj)
        variant_uid_to_count_dict[key] += 1

    # Now, get rid of rows that have no Sample associated, when there is more
    # than one row for that Variant.
    variant_uid_to_deleted_row_dict = {}
    modified_obj_list = []
    for obj in obj_list:
        key = _make_catch_all_key(obj)
        is_not_associated_with_sample = bool(
                not obj[MELTED_SCHEMA_KEY__ES_UID])
        exists_some_row_with_sample_association = bool(
                variant_uid_to_count_dict[key] > 1)
        if (is_not_associated_with_sample and
                exists_some_row_with_sample_association):
            assert key not in variant_uid_to_deleted_row_dict, (
                    "We've already seen the catch-all row for this Variant.")
            variant_uid_to_deleted_row_dict[key] = obj
        else:
            modified_obj_list.append(obj)

    # Now add back data for all the sets associated with each Variant.
    for obj in modified_obj_list:
        key = _make_catch_all_key(obj)
        if key in variant_uid_to_deleted_row_dict:
            catch_all_obj = variant_uid_to_deleted_row_dict[key]
            obj[ALL_VS_LABEL_KEY] = catch_all_obj[MELTED_SCHEMA_KEY__VS_LABEL]
            obj[ALL_VS_UID_KEY] = catch_all_obj[MELTED_SCHEMA_KEY__VS_UID]
        else:
            obj[ALL_VS_LABEL_KEY] = []
            obj[ALL_VS_UID_KEY] = []

    return modified_obj_list


def _prepare_additional_visible_keys(visible_key_list, key_to_parent_map):
    """Prepare all additional keys.
    """
    # Now we use this map to construct the field dictionary objects that can
    # be used for further processing.
    return [_prepare_visible_key_name_for_adapting_to_fe(key,
            key_to_parent_map) for key in visible_key_list]


def _prepare_visible_key_name_for_adapting_to_fe(key_name, key_to_parent_map):
    """Prepare single key for adapting to frontend.

    Returns:
        Dictionary representation of the field that can be handled by
        the next adaptation method like adapt_non_recursive(), e.g.:
        {
            'field': 'gt_nums',
            'is_subkey': True,
            'parent_col': 've_data'
        }
    """
    result = {'field': key_name}

    # Maybe add the parent col this field comes from, for properly handling
    # parsing.
    parent_col = key_to_parent_map.get(key_name, None)
    if parent_col is not None:
        result['is_subkey'] = True
        result['parent_col'] = parent_col

    # Hack: Prevent variant set uid from showing up in ui.
    if key_name == MELTED_SCHEMA_KEY__VS_UID:
        result['hide'] = True

    return result


def adapt_melted_object_field(val, fdict):
    """
    Any extra processing we have to do to the melted string according to
    fdict keys should happen here.
    """
    if 'recase' in fdict:
        val = titlecase_spaces(val)
    return val


def adapt_cast_object_list_field(cast_object_dict_list, fdict):
    """Converts a Cast object's field as a list into a string where values
    have been bucketed by unique type, with counts in parens.

    TODO: This doesn't really make sense for fields that take on continuous
    number values. Figure out what to do with these kinds of fields.
    """
    # First, extract the relevant fields and recase them.
    value_list = []
    for cast_object_dict in cast_object_dict_list:
        if cast_object_dict is None:
            continue
        else:
            val = cast_object_dict.get(fdict['field'], None)
            if val is None: continue
            if 'recase' in fdict:
                val = titlecase_spaces(val)
            value_list.append(val)

    # Using ordered dicts with None as an ordered set
    buckets = OrderedDict()
    for val in value_list:
        key = str(val)
        buckets[key] = buckets.get(key, 0) + 1

    # If gathering, just list all values, maintaining order.
    if fdict.get('format','bucket') is 'gather':
        # Create string.
        return (' | '.join(map(str, buckets.keys())))

    # If bucketing, count each.
    elif fdict.get('format','bucket') is 'bucket':
        # Create string.
        return (' | '.join(['%s (%d)' % (key, count)
                for key, count in buckets.iteritems()]))


def adapt_gene_list_to_frontend(obj_list, alignment_group):
    """
    Returns:
        JSON string with the objects in the form that can be drawn by the
        Datatables component.
    """
    # TODO: This is an initial implementation. It can be made prettier.
    gene_view_fields = ['gene', 'num_variants']

    obj_field_config = [{
        'mData': field.upper(),
        'sTitle': field
    } for field in gene_view_fields]


    fe_obj_list = []
    for obj in obj_list:
        gene_href = get_href_for_gene(alignment_group, obj['GENE'])
        fe_obj_list.append({
            'GENE': '<a href="' + gene_href + '">' + obj['GENE'] + '</a>',
            'NUM_VARIANTS': obj['NUM_VARIANTS']
        })

    return json.dumps({
        'obj_list': fe_obj_list,
        'field_config': obj_field_config
    })


def get_href_for_gene(alignment_group, gene):
    """Generates href to analyze variants view filtered by this gene.
    """
    analyze_tab_part = reverse(
            'main.views.tab_root_analyze',
            args=(alignment_group.reference_genome.project.uid,
                    alignment_group.uid,
                    'variants'))
    gene_filter_part = '?filter=INFO_EFF_GENE=' + gene
    return analyze_tab_part + gene_filter_part
