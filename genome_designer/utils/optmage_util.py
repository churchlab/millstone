"""Utils for connecting Millstone to optMAGE.
"""

from optmage.oligo_designer import OligoGenerator
from optmage.oligo_designer import OligoTarget
from optmage.oligo_designer import OligoWriter
from optmage.oligo_designer import OptMAGEConfig
from optmage.oligo_designer import OPT_MAGE_MUTATION_TYPE__SUBSTITUTION
from optmage.oligo_designer import OLIGO_TARGET_REQUIRED_PARAMS

from main.exceptions import ValidationException
from variants.vcf_parser import SV_TYPES


def print_mage_oligos(variant_set, output_filehandle, target_id_prefix):
    """Prints MAGE oligos for Variants in variant_set.

    Args:
        variant_set: A VariantSet.
        output_filehandle: Target for writing the designed oligos.

    Raises:
        ValidationException
    """
    # First validate that we support printing oligos for this type.
    _validate_variant_set_for_printing_mage_oligos(variant_set)

    config = OptMAGEConfig()

    oligo_target_list = []
    for variant in variant_set.variants.all():
        alt_value = variant.variantalternate_set.all()[0].alt_value
        oligo_target_list.append(OligoTarget(config, {
                'target_id': variant.position,
                'strand': 1,
                'start': variant.position,
                'end': variant.position + 1,
                'mutation_type': OPT_MAGE_MUTATION_TYPE__SUBSTITUTION,
                'mutation_seq': alt_value
        }))

    oligo_generator = OligoGenerator(config)

    oligo_result_list = [
            oligo_generator.generate_oligo(oligo_target)
            for oligo_target in oligo_target_list]
    OligoWriter.write_default(oligo_result_list, output_filehandle)


def _validate_variant_set_for_printing_mage_oligos(variant_set):
    """Validates the VariantSet to make sure we support printing oligos
    of this type.

    Raises:
        ValidationException
    """
    for variant in variant_set.variants.all():
        # Don't support SV types.
        if variant.type in SV_TYPES.values():
            raise ValidationException("SVs Not supported")

        # Don't support case of multiple alts.
        if variant.variantalternate_set.count() != 1:
            raise ValidationException("All Variants must have exactly one alt.")

        ref_value = variant.ref_value
        alt_value = variant.variantalternate_set.all()[0]

        if ref_value == alt_value:
            raise ValidationException("Nothing to do.")


# def _determine_mutation_type(variant):
#     if variant.ref ==
