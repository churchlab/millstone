"""
Module created, admittedly, out of failure to come up with better solution for
circular dependency tendency of design.
"""

from main.model_views import MeltedVariantView


def variant_as_melted_list(variant_obj, variant_id_to_metadata_dict=None):
    """Melt the variant into a list of objects, one per sample.

    Args:
        variant_obj: The Variant object to melt.
        variant_id_to_metadata_dict: Restrict melted entities to those
            related to these samples.

    Returns:
        List of MeltedVariantView objects corresponding to the variant.
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

    if variant_id_to_metadata_dict is not None:
        passing_sample_ids = variant_id_to_metadata_dict[variant_obj.id].get(
                'passing_sample_ids', set())
    else:
        passing_sample_ids = None

    # Iterate over the evidence objects.
    for variant_evidence_obj in common_data_obj.variantevidence_set.all():
        if passing_sample_ids is not None:
            sample_id = variant_evidence_obj.experiment_sample.id
            if not sample_id in passing_sample_ids:
                continue
        melted_list.append(MeltedVariantView(variant_obj, common_data_obj,
            variant_evidence_obj))

    return melted_list
