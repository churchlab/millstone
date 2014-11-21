"""Utility functions for testing sv callers.
"""

import os

from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from variants.vcf_parser import extract_raw_data_dict


def verify_variant_type(variants, variant_type, pos, length):
    """"Helper function for checking a specific variant type.

    Args:
        variants: List of variants to look through.
        variant_type: The type of variant to verify.
        pos: Position of the variant.
        length: Size of the variant.

    Returns:
        Boolean indicating whether variant was verified.
    """
    for variant in variants:
        # Check variant against following gauntlet.
        if variant['type'] != variant_type:
            continue # Fail, incorrect type.
        if abs(abs(variant['pos']) - pos) >= 100:
            continue # Fail, incorrect position.
        if (length != -1 and
                abs(abs(variant['length']) - length) >= 100):
            continue # Fail, incorrect length.
        # Success, variant made it through gauntlet.
        return True

    # If we got here, no matches were found.
    return False


def get_vcf_files(alignment_group):
    """Gets vcf files related to the AlignmentGroup.

    Returns:
        Dict mapping from vcf type to file location.
    """
    vcf_files = {}
    vcf_types = [VARIANT_TOOL_PARAMS_MAP[tool]['dataset_type']
            for tool in settings.ENABLED_VARIANT_CALLERS]
    for vcf_type in vcf_types:
        vcf_dataset = get_dataset_with_type(alignment_group, vcf_type)
        if vcf_dataset is None:
            continue
        vcf_location = vcf_dataset.get_absolute_location()
        assert os.path.exists(vcf_location)
        vcf_files[vcf_type] = vcf_location
    return vcf_files


def get_sv_variants(alignment_group, vcf_type):
    """Pull out SV variants.
    """
    vcf_files = get_vcf_files(alignment_group)
    assert vcf_type in vcf_files
    variants = []
    with open(vcf_files[vcf_type]) as fh:
        vcf_reader = vcf.Reader(fh)
        for record_idx, record in enumerate(vcf_reader):

            raw_data_dict = extract_raw_data_dict(record)

            # we should expect exactly 1 alternate
            assert len(raw_data_dict['INFO_SVLEN']) == 1, (
                'length of INFO_SVLEN > 1: {svlen}'.format(
                        svlen=raw_data_dict['INFO_SVLEN']))
            assert len(raw_data_dict['INFO_SVTYPE']) == 1, (
                'length of INFO_SVLEN > 1: {svtype}'.format(
                        svtype=raw_data_dict['INFO_SVTYPE']))

            variant_type = str(raw_data_dict.get('INFO_SVTYPE',
                    raw_data_dict.get('TYPE'))[0])
            pos = int(raw_data_dict.get('POS'))
            length = int(raw_data_dict.get('INFO_SVLEN')[0])
            variants.append({
                'type': variant_type,
                'pos': pos,
                'length': length
                })
    return variants
