"""
Utility methods for manipulating variant sets.
"""

from django.core.exceptions import ObjectDoesNotExist

from main.constants import UNDEFINED_STRING
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet

import re


MODIFY_VARIANT_SET_MEMBERSHIP__ADD = 'add'
MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE = 'remove'
VALID_ACTIONS = set([
    MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
    MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE
])

# Lenient regex that matches uids.
UID_REGEX = re.compile('\w+')


def update_variant_in_set_memberships(permissions_verified_ref_genome,
        uid_data_str_list, action, variant_set_uid):
    """Modifies the memberships of the given list of variant/sample objects
    in the given VariantSet.

    Args:
        permissions_verified_ref_genome: Used to verify ownership of contained
            entities. Clients of this method must confirm that the requestor
            has permissions to modify this ReferenceGenome.
        uid_data_str_list: A list of strings in a convoluted format that Gleb
            made up while hacking an initial solution. It is of the form:
                '<variant_uid>,<sample1_uid> | <sample2_uid>'
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
    ### Initial validation.

    # Make sure we have a valid action.
    if not action in VALID_ACTIONS:
        return {
            'alert_type': 'error',
            'alert_msg': 'Bad action type.'
        }

    # Make sure that the VariantSet is valid.
    try:
        variant_set = VariantSet.objects.get(
                reference_genome=permissions_verified_ref_genome,
                uid=variant_set_uid)
    except ObjectDoesNotExist as e:
        return {
            'alert_type': 'error',
            'alert_msg': 'Variant Set does not exist or insufficient permissions.'
        }

    ### Parse the convoluted format.

    grouped_uid_dict_list = []
    for data_str in uid_data_str_list:
        variant_uid_and_maybe_sample_uid_list = data_str.split(',')
        variant_uid = variant_uid_and_maybe_sample_uid_list[0]
        sample_uid_list = []
        if len(variant_uid_and_maybe_sample_uid_list) > 1:
            assert len(variant_uid_and_maybe_sample_uid_list) == 2, (
                    "Actual: %d" % len(variant_uid_and_maybe_sample_uid_list))
            sample_uid_list = variant_uid_and_maybe_sample_uid_list[1].split('|')
            sample_uid_list = [raw_sample_uid.strip() for raw_sample_uid
                    in sample_uid_list]
            for sample_uid in sample_uid_list:
                if UID_REGEX.match(sample_uid):
                    grouped_uid_dict_list.append({
                        'variant_uid': variant_uid,
                        'sample_uid': sample_uid
                    })
                elif sample_uid == UNDEFINED_STRING:
                    grouped_uid_dict_list.append({
                        'variant_uid': variant_uid,
                        'sample_uid': UNDEFINED_STRING
                    })
                else:
                    raise AssertionError(
                            "Unknown sample_uid format: %s" % sample_uid)

    ### Query the relevant objects.

    unique_variant_uid_set = set(
            [group['variant_uid'] for group in grouped_uid_dict_list])
    variant_uid_to_obj_map = _get_uid_to_db_object_map(Variant,
            unique_variant_uid_set,
            {'reference_genome': permissions_verified_ref_genome})

    unique_sample_uid_set = set([group['sample_uid'] for group
            in grouped_uid_dict_list
            if not group['sample_uid'] == UNDEFINED_STRING])
    sample_uid_to_obj_map = _get_uid_to_db_object_map(ExperimentSample,
            unique_sample_uid_set,
            {'project': permissions_verified_ref_genome.project})

    # Perform the appropriate action.
    if action == MODIFY_VARIANT_SET_MEMBERSHIP__ADD:
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
    else:
        for group in grouped_uid_dict_list:
            variant = variant_uid_to_obj_map[group['variant_uid']]
            try:
                vtvs = VariantToVariantSet.objects.get(
                        variant=variant,
                        variant_set=variant_set)
            except ObjectDoesNotExist as e:
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
            except ObjectDoesNotExist as e:
                # It's possible this relation was removed on an earlier
                # iteration.
                # TODO(gleb): Think about what could go wrong here and add
                # tests.
                continue
            sample_uid = group['sample_uid']
            if sample_uid == UNDEFINED_STRING:
                if vtvs.sample_variant_set_association.count() == 0:
                    vtvs.delete()

    # Return success response if we got here.
    return {
        'alert_type': 'info',
        'alert_msg': 'success'
    }


def _get_uid_to_db_object_map(model_class, unique_uid_set,
        query_kwargs={}):
    """Queries the database for the given list of uids, returning a map from
    uid to object.
    """
    uid_to_obj_map = {}
    query_set = model_class.objects.filter(**query_kwargs).filter(
            uid__in=unique_uid_set)
    for obj in query_set:
        uid_to_obj_map[obj.uid] = obj
    return uid_to_obj_map
