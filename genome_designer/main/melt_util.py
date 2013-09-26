"""
Module created, admittedly, out of failure to come up with better solution for
circular dependency tendency of design.
"""

from main.model_views import MeltedVariantView


def variant_as_melted_list(variant_obj, variant_id_to_metadata_dict=None):
    """Melt the variant into a list of objects, one per common data per sample.

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
            passing_sample_ids = variant_id_to_metadata_dict[variant_obj.id].get(
                    'passing_sample_ids', set())
        else:
            passing_sample_ids = None

        for variant_evidence_obj in variant_evidence_list:
            if passing_sample_ids is not None:
                sample_id = variant_evidence_obj.experiment_sample.id
                if not sample_id in passing_sample_ids:
                    continue
            melted_list.append(MeltedVariantView(variant_obj, common_data_obj,
                    variant_evidence_obj))

    return melted_list
