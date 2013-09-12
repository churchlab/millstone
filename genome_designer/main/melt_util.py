"""
Module created, admittedly, out of failure to come up with better solution for
circular dependency tendency of design.
"""

from main.model_views import MeltedVariantView


def variant_as_melted_list(variant_obj):
    """Returns a list of MeltedVariantView objects corresponding to the variant.
    """
    melted_list = []

    # Get the common data object.
    all_common_data = variant_obj.variantcallercommondata_set.all()
    if len(all_common_data) == 0:
        # TODO: Confirm this default behavior is correct.
        return [variant_obj]
    if len(all_common_data) > 1:
        raise AssertionError("Support for multiple VariantCallerCommonData " +
                "objects not implemented yet.")
    common_data_obj = all_common_data[0]

    if len(common_data_obj.variantevidence_set.all()) == 0:
        return [variant_obj]

    # Iterate over the evidence objects.
    for variant_evidence_obj in common_data_obj.variantevidence_set.all():
        melted_list.append(MeltedVariantView(variant_obj, common_data_obj,
            variant_evidence_obj))

    return melted_list
