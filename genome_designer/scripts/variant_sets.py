from django.http import HttpResponseBadRequest

from main.models import AlignmentGroup
from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet

def add_or_remove_variants_from_set(variantUidList, variantSetAction,
        variantSetUid, sampleUidList=[]):

    # Parse the data and look up the relevant model instances:

    # Variant list.
    variant_list = Variant.objects.filter(
            uid__in=variantUidList)
    assert len(variant_list) == len(variantUidList)
    if not len(variantSetUid) > 0:
        response_json = {
                'alert_type':'error',
                'alert_msg':'At least one variant required.'
        }
        return response_json

    # Variant set.
    try:
        variant_set = VariantSet.objects.get(uid=variantSetUid)
    except ObjectDoesNotExist as e:
        response_json = {
                'alert_type' : 'error',
                'alert_msg' : 'Variant set does not exist.'
        }
        return response_json

    # Identify existing variant to set relationships.
    vvs_existing = VariantToVariantSet.objects.filter(
            variant__in=variant_list,
            variant_set=variant_set)

    exist_count = len(vvs_existing)
    total_count = len(variant_list)

    # Perform the variant set action.
    if variantSetAction == 'add':


        # Identify existing relationships if any exist.
        if len(vvs_existing) > 0:

            var_already_in_set = variant_list.filter(
                    variantset=variant_set)

            response_json = {
                    'alert_type' : 'warn',
                    'alert_msg' : ('%d of %d variants were already in the' +
                            ' chosen set and were ignored.') % (exist_count,
                            total_count)
            }

        # No existing relationships.
        else:
            var_already_in_set = []
            response_json = {
                    'alert_type' : 'info',
                    'alert_msg' : '%d variants successfully added.' % (
                            total_count)
            }

        # Create new vvs objects, subtracting existing.
        vvs_to_create = [VariantToVariantSet(
                variant=var,
                variant_set=variant_set)
                for var in variant_list if var not in var_already_in_set]

        VariantToVariantSet.objects.bulk_create(vvs_to_create)
        return response_json

    # Perform the variant remove action.
    elif variantSetAction == 'remove':

        if len(vvs_existing) == 0:
            response_json = {
                    'alert_type' : 'error',
                    'alert_msg' : 'None of the selected variants are in the' +
                            ' chosen set.'
            }
            return response_json

        elif len(vvs_existing) != len(variant_list):
            vvs_existing.delete()

            extra_count = total_count - exist_count

            response_json = {
                    'alert_type':'warn',
                    'alert_msg': ('%d of %d variants were not in the chosen' +
                            ' set and were ignored.') % (extra_count,
                            total_count)
            }
            return response_json

        else:
            vvs_existing.delete()
            response_json = {
                    'alert_type':'info',
                    'alert_msg':'%d variants successfully removed.' % (
                            total_count)
            }
            return response_json
    else:
        response_json = {
                'alert_type' : 'error',
                'alert_msg': 'Bad variantSet action type.'
        }
        return response_json
