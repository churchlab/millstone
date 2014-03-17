"""
Utility methods for manipulating variant sets.
"""

import re

from django.core.exceptions import ObjectDoesNotExist

from main.constants import UNDEFINED_STRING
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet


MODIFY_VARIANT_SET_MEMBERSHIP__ADD = 'add'
MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE = 'remove'
VALID_ACTIONS = set([
    MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
    MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE
])

# Lenient regex that matches uids.
UID_REGEX = re.compile('\w+')


def update_variant_in_set_memberships(ref_genome, uid_data_str_list,
        action, variant_set_uid):
    """Modifies the memberships of the given list of variant/sample objects
    in the given VariantSet.

    Args:
        ref_genome: Used to verify ownership of contained
            entities. Clients of this method must confirm that the requestor
            has permissions to modify this ReferenceGenome.
        uid_data_str_list: A list of strings, each of one of two possible forms:
            * "variant_uid"
            * "variant_uid, sample_uid"
        action: The action to perform.
        variant_set_uid: The set to the add the variant to.

    Expected exceptions are caught and reported in the response object.

    Returns:
        A dictionary containing information about how the request was handled.
        This is ultimately returned to the UI to show the client a message.
        Contains the following keys:
            * alert_type: Type of message. Either 'info', 'error', or 'warn'.
            * alert_msg: Additional information shown to the user.
    """
    validation_result = _initialial_validation(
            uid_data_str_list, action)
    if validation_result['alert_type'] == 'error':
        return validation_result

    # Make sure that the VariantSet is valid.
    try:
        variant_set = VariantSet.objects.get(
                reference_genome=ref_genome,
                uid=variant_set_uid)
    except ObjectDoesNotExist:
        return {
            'alert_type': 'error',
            'alert_msg': 'Variant Set does not exist or insufficient permissions.'
        }

    # Convert to object interface to help clarify further processing.
    grouped_uid_dict_list = _convert_pair_list_into_object_list(
            uid_data_str_list)

    # Get helper objects for the query.
    (variant_uid_to_obj_map, sample_uid_to_obj_map) = (
            _get_cached_uid_to_object_maps(ref_genome, grouped_uid_dict_list))

    # Perform modification.
    if action == MODIFY_VARIANT_SET_MEMBERSHIP__ADD:
        _perform_add(grouped_uid_dict_list, variant_set,
                variant_uid_to_obj_map, sample_uid_to_obj_map)
    else: # action == MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE
        _perform_remove(grouped_uid_dict_list, variant_set,
                variant_uid_to_obj_map, sample_uid_to_obj_map)

    # These actions invalidate the materialized view.
    ref_genome.invalidate_materialized_view()

    # Return success response if we got here.
    return {
        'alert_type': 'info',
        'alert_msg': 'success'
    }


def _initialial_validation(uid_data_str_list, action):
    """Initial validation, or why statically compiled languages have their
    upsides.
    """
    # Make sure we have a valid action.
    if not action in VALID_ACTIONS:
        return {
            'alert_type': 'error',
            'alert_msg': 'Bad action type: %s' % action
        }

    # Make sure uid_data_str_list is properly formatted.
    for uid_data_str in uid_data_str_list:
        parts = uid_data_str.split(',')
        if not (0 < len(parts) <= 2):
            return {
                'alert_type': 'error',
                'alert_msg': 'Bad variant uid / sample uid pair: %s' % uid_data_str
            }

    # All good. Return filler message.
    return {
        'alert_type': 'info',
        'alert_msg': 'Validation passed.'
    }


def _convert_pair_list_into_object_list(uid_data_str_list):
    """Converts the list of pairs into a list of objects with keys:
    * variant_uid
    * sample_uid
    """
    grouped_uid_dict_list = []
    for uid_data_str in uid_data_str_list:
        parts = uid_data_str.split(',')
        variant_uid = parts[0]
        if len(parts) > 1:
            sample_uid = parts[1]
        else:
            sample_uid = ''
        if UID_REGEX.match(sample_uid):
            grouped_uid_dict_list.append({
                'variant_uid': variant_uid,
                'sample_uid': sample_uid
            })
        else:
            grouped_uid_dict_list.append({
                'variant_uid': variant_uid,
                'sample_uid': UNDEFINED_STRING
            })
    return grouped_uid_dict_list


def _get_uid_to_db_object_map(model_class, unique_uid_set,
        query_kwargs={}):
    """Queries the database for the given list of uids, returning a map from
    uid to Python object.
    """
    uid_to_obj_map = {}
    query_set = model_class.objects.filter(**query_kwargs).filter(
            uid__in=unique_uid_set)
    for obj in query_set:
        uid_to_obj_map[obj.uid] = obj
    return uid_to_obj_map


def _get_cached_uid_to_object_maps(ref_genome, grouped_uid_dict_list):
    """Reduces the number of DB reads necessary to get the associated entities.
    """
    unique_variant_uid_set = set(
            [group['variant_uid'] for group in grouped_uid_dict_list])
    variant_uid_to_obj_map = _get_uid_to_db_object_map(Variant,
            unique_variant_uid_set,
            {'reference_genome': ref_genome})

    unique_sample_uid_set = set([group['sample_uid'] for group
            in grouped_uid_dict_list
            if not group['sample_uid'] == UNDEFINED_STRING])
    sample_uid_to_obj_map = _get_uid_to_db_object_map(ExperimentSample,
            unique_sample_uid_set,
            {'project': ref_genome.project})
    return (variant_uid_to_obj_map, sample_uid_to_obj_map)


def _perform_add(grouped_uid_dict_list, variant_set, variant_uid_to_obj_map,
        sample_uid_to_obj_map):
    for group in grouped_uid_dict_list:
        variant = variant_uid_to_obj_map[group['variant_uid']]
        vtvs, created = VariantToVariantSet.objects.get_or_create(
                variant=variant,
                variant_set=variant_set)
        # Maybe add sample association.
        sample_uid = group['sample_uid']
        if sample_uid == UNDEFINED_STRING:
            continue
        vtvs.sample_variant_set_association.add(
                sample_uid_to_obj_map[sample_uid])


def _perform_remove(grouped_uid_dict_list, variant_set,
        variant_uid_to_obj_map, sample_uid_to_obj_map):
    for group in grouped_uid_dict_list:
        variant = variant_uid_to_obj_map[group['variant_uid']]
        try:
            vtvs = VariantToVariantSet.objects.get(
                    variant=variant,
                    variant_set=variant_set)
        except ObjectDoesNotExist:
            # It's possible this relation was removed on an earlier
            # iteration.
            # TODO(gleb): Think about what could go wrong here and add
            # tests.
            continue
        sample_uid = group['sample_uid']
        if sample_uid == UNDEFINED_STRING:
            # See next for loop.
            continue
        vtvs.sample_variant_set_association.remove(
                sample_uid_to_obj_map[sample_uid])

    # Iterate again and destroy any VariantToVariantSet that
    # should be completely removed.
    for group in grouped_uid_dict_list:
        variant = variant_uid_to_obj_map[group['variant_uid']]
        try:
            vtvs = VariantToVariantSet.objects.get(
                    variant=variant,
                    variant_set=variant_set)
        except ObjectDoesNotExist:
            # It's possible this relation was removed on an earlier
            # iteration.
            # TODO(gleb): Think about what could go wrong here and add
            # tests.
            continue
        sample_uid = group['sample_uid']
        if sample_uid == UNDEFINED_STRING:
            if vtvs.sample_variant_set_association.count() == 0:
                vtvs.delete()
