"""Utils for connecting Millstone to optMAGE.
"""

from optmage.oligo_designer import OligoGenerator
from optmage.oligo_designer import OligoTarget
from optmage.oligo_designer import OligoWriter
from optmage.oligo_designer import OptMAGEConfig
from optmage.oligo_designer import OPT_MAGE_MUTATION_TYPE__REFERENCE
from optmage.oligo_designer import OPT_MAGE_MUTATION_TYPE__SUBSTITUTION

from main.exceptions import ValidationException
from variants.vcf_parser import SV_TYPES


# Whether the oligo introduces or reverts mutation.
OPT_MAGE_EXPERIMENT_DIR__FORWARD = 'F'
OPT_MAGE_EXPERIMENT_DIR__REVERSE = 'R'
OPT_MAGE_VALID_EXPERIMENT_DIRS = (
    OPT_MAGE_EXPERIMENT_DIR__FORWARD,
    OPT_MAGE_EXPERIMENT_DIR__REVERSE
)


def print_mage_oligos(variant_set, output_filehandle, target_id_prefix,
        experiment_dir=OPT_MAGE_EXPERIMENT_DIR__FORWARD):
    """Prints MAGE oligos for Variants in variant_set.

    Args:
        variant_set: A VariantSet.
        output_filehandle: Target for writing the designed oligos.
        target_id_prefix: String prefix used for target id string in output.
        experiment_dir: Direction to print oligos. Choose from VALID_PRINT_OLIGOS_DIR.

    Raises:
        ValidationException
    """
    assert experiment_dir in OPT_MAGE_VALID_EXPERIMENT_DIRS

    # First validate that we support printing oligos for this type.
    passing_variants = _validate_variant_set_for_printing_mage_oligos(
            variant_set)

    config = OptMAGEConfig()

    oligo_target_list = []
    for variant in passing_variants:
        alt_value = variant.variantalternate_set.all()[0].alt_value

        # Configure the oligo for this variant.
        target_params = {
            'target_id': variant.position,
            'strand': 1,
            'start': variant.position,
            'end': variant.position + 1,
        }
        if experiment_dir == OPT_MAGE_EXPERIMENT_DIR__REVERSE:
            target_params.update({
                'mutation_type': OPT_MAGE_MUTATION_TYPE__REFERENCE,
            })
        else:
            target_params.update({
                'mutation_type': OPT_MAGE_MUTATION_TYPE__SUBSTITUTION,
                'mutation_seq': alt_value,
            })
        oligo_target_list.append(OligoTarget(config, target_params))

    oligo_generator = OligoGenerator(config)

    oligo_result_list = [
            oligo_generator.generate_oligo(oligo_target)
            for oligo_target in oligo_target_list]
    OligoWriter.write_default(oligo_result_list, output_filehandle)


def _validate_variant_set_for_printing_mage_oligos(variant_set):
    """Validates the VariantSet to make sure we support printing oligos
    of this type. Discards Variants for which it doesn't make sense to
    make oligos (e.g if no VariantAlternate).

    Args:
        variant_set: The VariantSet we are printing oligos for.

    Returns:
        List of Variants that pass.

    Raises:
        ValidationException
    """
    passing_variants = []
    for variant in variant_set.variants.all():
        # Don't support SV types.
        if variant.type in SV_TYPES.values():
            raise ValidationException("SVs Not supported")

        # No alt, silently skip.
        if variant.variantalternate_set.count() == 0:
            continue

        # Don't support case of multiple alts.
        if variant.variantalternate_set.count() > 1:
            raise ValidationException(
                    "All Variants must have exactly one alt. " +
                    "Variant with uid " + variant.uid + " has " +
                    str(variant.variantalternate_set.count()))

        ref_value = variant.ref_value
        alt_value = variant.variantalternate_set.all()[0]

        if ref_value == alt_value:
            raise ValidationException("Nothing to do.")

        passing_variants.append(variant)
    return passing_variants
