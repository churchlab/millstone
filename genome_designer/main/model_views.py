"""
Classes that describe how a particular model should be viewed.
"""

import json
import string

from django.core.urlresolvers import reverse
from django.db.models.query import QuerySet

from main.constants import UNDEFINED_STRING
from main.data_util import cast_object_list_field_as_bucket_string
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantAlternate
from main.models import VariantEvidence
from main.models import VariantSet
from scripts.dynamic_snp_filter_key_map import MAP_KEY__COMMON_DATA
from scripts.dynamic_snp_filter_key_map import MAP_KEY__ALTERNATE
from scripts.dynamic_snp_filter_key_map import MAP_KEY__EVIDENCE
from variants.common import VARIANT_KEY_TO_MATERIALIZED_VIEW_COL_MAP


class BaseVariantView(object):
    """Common methods for model views.
    """

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        if 'visible_key_names' in kwargs:
            variant_key_map = kwargs['variant_key_map']
            visible_key_names = kwargs['visible_key_names']

            additional_common_data_field_list = [field for field in
                    visible_key_names if field in
                    variant_key_map[MAP_KEY__COMMON_DATA]]

            additional_alternate_field_list = [field for field in
                    visible_key_names if field in
                    variant_key_map[MAP_KEY__ALTERNATE]]

            additional_evidence_field_list = [field for field in visible_key_names
                    if field in variant_key_map[MAP_KEY__EVIDENCE]]
        else:
            additional_common_data_field_list = []
            additional_evidence_field_list = []
            additional_alternate_field_list = []

        # Sort the visible keys into appropriate buckets based on the map.
        return (Variant.get_field_order() +
                VariantAlternate.get_field_order(
                        additional_field_list=additional_alternate_field_list) + 
                VariantCallerCommonData.get_field_order(
                        additional_field_list=additional_common_data_field_list) +
                VariantEvidence.get_field_order(
                        additional_field_list=additional_evidence_field_list))


class CastVariantView(BaseVariantView):
    """View of Variant that does not expose the individual
    VariantToExperimentSample relationships.
    """
    def __init__(self, variant, data_cache):
        """Constructor.

        Args:
            variant: A Variant instance.
            data_cache: RequestScopedVariantDataCache instance.
        """
        self.variant = variant
        self.variantalternate_list = (
                data_cache.get_variant_alternate_list(variant))
        self.variant_caller_common_data_list = (
                data_cache.get_variant_caller_common_data_list(variant))
        self.variant_evidence_list = (
                data_cache.get_variant_evidence_list(variant))

    def custom_getattr(self, attr):
        """For an attribute of a Variant, returns what you would expect.

        For many-to-one relations, e.g. VariantCallerCommonData and
        VariantEvidence, returns a '|'-separated list of values for
        that attribute.

        Returns:
            A string representing the view for the attribute.
        """
        delegate_order = [
                self.variant,
                self.variantalternate_list,
                self.variant_caller_common_data_list,
                self.variant_evidence_list
        ]
        for delegate in delegate_order:
            result_list = []

            # Convert to list for next step.
            if isinstance(delegate, QuerySet) or isinstance(delegate, list):
                delegate_list = delegate
            else:
                delegate_list = [delegate]

            # Iterate through the delegates.
            for delegate in delegate_list:
                if hasattr(delegate, attr):
                    result_list.append(getattr(delegate, attr))
                else:
                    try: result_list.append(delegate[attr])
                    except: pass

            # If no results, continue to next delegate.
            if not len(result_list):
                continue

            # If exactly one object, just return that object.
            if len(result_list) == 1:
                return result_list[0]

            # If 4 or less, return a '|'-separated string.
            # Number 4 is about how many short strings can reasonably fit in a
            # cell in the ui.
            # NOTE: This case is most useful for showing Variant ALTs.
            if len(result_list) <= 4:
                return ' | '.join([str(res) for res in result_list])

            # Otherwise return the count.
            return '{' + str(len(result_list)) + '}'

        # Default.
        return UNDEFINED_STRING

    @classmethod
    def variant_as_cast_view(clazz, variant_obj, variant_id_to_metadata_dict,
            data_cache):
        """Factory method returns a cast view for the given variant.

        Args:
            variant_obj: The Variant object to melt.
            variant_id_to_metadata_dict: See TODO.
            data_cache: Object that handles bulk SQL lookup.

        Returns:
            A CastVariantView instance.
        """
        # TODO: Do we need to modify the view based on passing samples in the
        # metadata.
        return CastVariantView(variant_obj, data_cache)


class MeltedVariantView(BaseVariantView):
    """View of a Variant when showing one row per VariantToExperimentSample.
    """

    def __init__(self, variant, variant_caller_common_data=None,
            variant_evidence=None):
        self.variant = variant
        self.variantalternate_list = VariantAlternate.objects.filter(
                variant=variant, variantevidence=variant_evidence)
        self.variant_caller_common_data = variant_caller_common_data
        self.variant_evidence = variant_evidence

    def custom_getattr(self, attr):
        """Custom implementation of getting an attribute.

        We need this since to make sure we delegate to the objects that compose
        this object, in the correct order.
        """
        # Handle certain fields manually.
        # TODO: Come up with more generic way of manually handling fields.
        if attr == 'variantset_set':
            # TODO(gleb): Make the queries below efficient. Right now I'm just
            # trying to get it to work.
            variant_set_list = []
            if self.variant_evidence is None:
                # Return VariantSets which contain no ExperimentSample
                # association.
                no_sample_variant_sets = []
                for vtvs in self.variant.varianttovariantset_set.all():
                    if len(vtvs.sample_variant_set_association.all()) == 0:
                        no_sample_variant_sets.append(vtvs.variant_set)
                variant_set_list = no_sample_variant_sets

            else:
                # Only show sets associated with the VariantEvidence.
                sample = self.variant_evidence.experiment_sample
                variant_sets_for_this_sample = []
                for vtvs in self.variant.varianttovariantset_set.all():
                    if sample in vtvs.sample_variant_set_association.all():
                        variant_sets_for_this_sample.append(vtvs.variant_set)
                variant_set_list = variant_sets_for_this_sample
            # TODO: Return a QuerySet more efficiently. This is a temporary bug fix.
            variant_set_uid_list = [vs.uid for vs in variant_set_list]
            return VariantSet.objects.filter(uid__in=variant_set_uid_list)

        delegate_order = [
                self.variant,
                self.variantalternate_list,
                self.variant_caller_common_data,
                self.variant_evidence
        ]

        variant_caller_common_data_map = (
                self.variant.reference_genome.get_variant_caller_common_map())

        for delegate in delegate_order:
            result_list = []

            # Convert to list for next step.
            if isinstance(delegate, QuerySet) or isinstance(delegate, list):
                delegate_list = delegate
            else:
                delegate_list = [delegate]

            # Iterate through the delegates.
            for delegate in delegate_list:

                if hasattr(delegate, attr):
                    result_list.append(getattr(delegate, attr))
                else:
                    try:
                        result_list.append(delegate[attr])
                    except:
                        pass

            if len(result_list) == 1:
                # If one object, return it.
                return result_list[0]
            elif len(result_list) > 1:
                # Else if more than one, return a '|'-separated string.
                return ' | '.join([str(res) for res in result_list])

        return UNDEFINED_STRING

    @classmethod
    def variant_as_melted_list(clazz, variant_obj,
            variant_id_to_metadata_dict=None):
        """Factory method that melts the variant into a list of objects,
        one per common data per sample.

        Args:
            variant_obj: The Variant object to melt.
            variant_id_to_metadata_dict: Restrict melted entities to those
                related to these samples.

        Returns:
            List of MeltedVariantView objects corresponding to the Variant.
            This list will always be at least length 1.
        """
        melted_list = []

        all_common_data = variant_obj.variantcallercommondata_set.all()

        # We guarantee returning at least one row.
        if len(all_common_data) == 0:
            melted_list.append(MeltedVariantView(variant_obj, None, None))
            return melted_list

        # Otherwise iterate through
        for common_data_obj in all_common_data:
            # Return one row not associated with any sample.
            melted_list.append(MeltedVariantView(variant_obj, common_data_obj, None))

            variant_evidence_list = common_data_obj.variantevidence_set.all()
            if len(variant_evidence_list) == 0:
                continue

            if variant_id_to_metadata_dict is not None:
                passing_sample_ids = (
                        variant_id_to_metadata_dict[variant_obj.id].get(
                                'passing_sample_ids', set()))
            else:
                passing_sample_ids = None

            for variant_evidence_obj in variant_evidence_list:
                if passing_sample_ids is not None:
                    # Don't add a row for VariantEvidence object if sample is
                    # not present.
                    sample_id = variant_evidence_obj.experiment_sample_id
                    if not sample_id in passing_sample_ids:
                        continue
                melted_list.append(MeltedVariantView(variant_obj, common_data_obj,
                        variant_evidence_obj))

        return melted_list


class GeneView(object):
    """View of a Gene.
    """

    def __init__(self, region):
        self.region = region
        self.label = region.label
        self.uid = region.uid

        # Assume that GENE region is composed of single interval.
        gene_interval = region.regioninterval_set.all()[0]
        self.start = gene_interval.start
        self.end = gene_interval.end

        # Count of Variants that occur in region.
        self.variants = Variant.objects.filter(
                reference_genome=region.reference_genome,
                position__gte=self.start,
                position__lt=self.end).count()

    @property
    def href(self):
        analyze_tab_part = reverse(
                'main.views.tab_root_analyze',
                args=(self.region.reference_genome.project.uid,
                        self.region.reference_genome.uid,
                        'variants'))
        gene_filter_part = '?filter=IN_GENE(' + self.label + ')'
        return analyze_tab_part + gene_filter_part


    @classmethod
    def get_field_order(clazz, **kwargs):
        return [
            {'field':'uid'},
            {'field':'label'},
            {'field':'start'},
            {'field':'end'},
            {'field':'variants'},
        ]


def adapt_non_recursive(obj_list, field_dict_list, reference_genome=None):
    """Adapts of list of objects that doesn't require recursive calling.

    Returns:
        JSON string representation of frontend objects.
    """
    # Parse the list of field names.
    field_list = [fdict['field'] for fdict in field_dict_list]

    # Get the verbose names for the fields.
    def _get_verbose(fdict):
        if 'verbose' in fdict:
            return fdict['verbose']
        else:
            return string.capwords(fdict['field'], '_').replace('_', ' ')
    field_verbose_names = [_get_verbose(fdict) for fdict in field_dict_list]

    # Aggregate list of objects that are ready for display by the frontend.
    fe_obj_list = []
    for melted_variant_obj in obj_list:
        visible_field_pairs = []
        for fdict in field_dict_list:
            value = None
            field = fdict['field']

            # HACK: Special handling for certain fields.
            if field == 'label':
                value = _create_single_variant_page_link_for_variant_object(
                        melted_variant_obj, reference_genome)
            elif field == 'position':
                value = _create_jbrowse_link_for_variant_object(
                        melted_variant_obj, reference_genome)

            elif fdict.get('is_subkey', False):
                assert 'parent_col' in fdict
                parent_dict_or_list = melted_variant_obj.get(fdict['parent_col'], {})
                if isinstance(parent_dict_or_list, dict):
                    # melted
                    parent_dict = parent_dict_or_list
                    value = str(parent_dict.get(field, ''))
                elif isinstance(parent_dict_or_list, list):
                    # cast
                    value = cast_object_list_field_as_bucket_string(
                            parent_dict_or_list, field)
            else:
                value = melted_variant_obj.get(field, '')

            # Pass empty string when no value present.
            if value is None:
                value = ''

            # Append (key, value) to list so we can create dict at the end.
            # Order doesn't matter here since it's determined in the
            # field config object constructed separately below.
            visible_field_pairs.append((field, value))

        fe_obj_list.append(dict(visible_field_pairs))

    # Create the config dict required by DataTables.js.
    obj_field_config = [{
        'mData': name,
        'sTitle': verbose_name
    } for (name, verbose_name) in zip(field_list, field_verbose_names)]

    return json.dumps({
        'obj_list': fe_obj_list,
        'field_config': obj_field_config
    })


def _create_single_variant_page_link_for_variant_object(variant_as_dict,
        reference_genome):
    """Constructs the label as an anchor that links to the single variant view.
    """
    full_href = reverse('main.views.single_variant_view',
                args=(reference_genome.project.uid, reference_genome.uid,
                        variant_as_dict['uid']))
    # Generate label from variant data.
    label = 'mut_' + str(variant_as_dict['position'])
    return '<a href="' + full_href + '">' + label + '</a>'


def _create_jbrowse_link_for_variant_object(variant_as_dict, reference_genome):
    """Constructs a JBrowse link for the Variant.
    """
    assert 'position' in variant_as_dict
    position = variant_as_dict['position']
    ref_genome_jbrowse_link = reference_genome.get_client_jbrowse_link()
    location_param = '&loc=' + str(position)
    full_href = ref_genome_jbrowse_link + location_param
    # TODO: Add alignment track param if this is a sample-specific view.
    return ('<a href="' + full_href + '" target="_blank">' + str(position) +
            '</a>')


MELTED_VARIANT_FIELD_DICT_LIST = [
    {'field': 'label', 'verbose': 'label'},
    {'field': 'position'},
    {'field': 'ref'},
    {'field': 'alt'},
    {'field': 'variant_set_label'},
    {'field': 'experiment_sample_uid'},
]

CAST_VARIANT_FIELD_DICT_LIST = [
    {'field': 'label', 'verbose': 'label'},
    {'field': 'position'},
    {'field': 'ref'},
    {'field': 'alt'},
    {'field': 'variant_sets'},
    {'field': 'total_samples'},
]


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
        default_field_dict_list = MELTED_VARIANT_FIELD_DICT_LIST
    else:
        default_field_dict_list = CAST_VARIANT_FIELD_DICT_LIST

    # Prepare additional fields for adaptation.
    additional_visible_field_dict_list = _prepare_additional_visible_keys(
            visible_key_names, reference_genome)

    all_field_dict_list = (default_field_dict_list +
            additional_visible_field_dict_list)

    return adapt_non_recursive(obj_list, all_field_dict_list, reference_genome)


def _prepare_additional_visible_keys(visible_key_list, reference_genome):
    """Prepare all additional keys.
    """
    # Some of the keys are actually inside of catch-all data objects. Our
    # materialized view has separate columns for each of these, so we need
    # a structure that maps from key name to column name.
    # NOTE: We enforce unique keys when dynamically generating
    # ReferenceGenome.variant_key_map.
    key_to_parent_map = _generate_key_to_materialized_view_parent_col(
            reference_genome)

    # Now we use this map to construct the field dictionary objects that can
    # be used for further processing.
    return [_prepare_visible_key_name_for_adapting_to_fe(key,
            key_to_parent_map) for key in visible_key_list]


def _generate_key_to_materialized_view_parent_col(reference_genome):
    key_to_parent_map = {}
    for submap_name, submap in reference_genome.variant_key_map.iteritems():
        for key in submap.iterkeys():
            assert not key in key_to_parent_map
            key_to_parent_map[key] = (
                    VARIANT_KEY_TO_MATERIALIZED_VIEW_COL_MAP.get(
                            submap_name, None))
    return key_to_parent_map


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
    result = {
        'field': key_name,
    }
    parent_col = key_to_parent_map.get(key_name, None)
    if parent_col is not None:
        result['is_subkey'] = True
        result['parent_col'] = parent_col
    return result
