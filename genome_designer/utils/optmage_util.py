"""Utils for connecting Millstone to optMAGE.
"""

from django.conf import settings
from optmage.oligo_designer import DEFAULT_REPLICATION_ORIGIN
from optmage.oligo_designer import DEFAULT_REPLICATION_TERMINUS
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


class ReplicationOriginParams(object):
    """Param object passed to print_mage_oligos.
    """

    def __init__(self, origin_start, origin_end, terminus_start, terminus_end):
        self.origin_start = int(origin_start)
        self.origin_end = int(origin_end)
        self.terminus_start = int(terminus_start)
        self.terminus_end = int(terminus_end)

    @classmethod
    def from_defaults(cls):
        """Factory method to create params from default.
        """
        return cls(
                DEFAULT_REPLICATION_ORIGIN[0],
                DEFAULT_REPLICATION_ORIGIN[1],
                DEFAULT_REPLICATION_TERMINUS[0],
                DEFAULT_REPLICATION_TERMINUS[1])

    def get_origin_interval(self):
        """Returns tuple pair.
        """
        return (self.origin_start, self.origin_end)

    def get_terminus_interval(self):
        """Returns tuple pair.
        """
        return (self.terminus_start, self.terminus_end)


def print_mage_oligos(variant_set, output_filehandle, target_id_prefix,
        replication_origin_params,
        experiment_dir=OPT_MAGE_EXPERIMENT_DIR__FORWARD):
    """Prints MAGE oligos for Variants in variant_set.

    Args:
        variant_set: A VariantSet.
        output_filehandle: Target for writing the designed oligos.
        target_id_prefix: String prefix used for target id string in output.
        replication_origin_params: ReplicationOriginParams object.
        experiment_dir: Direction to print oligos. Choose from VALID_PRINT_OLIGOS_DIR.

    Raises:
        ValidationException
    """
    assert experiment_dir in OPT_MAGE_VALID_EXPERIMENT_DIRS

    # First validate that we support printing oligos for this type.
    passing_variants = _validate_variant_set_for_printing_mage_oligos(
            variant_set)

    config = OptMAGEConfig()
    config.replication_origin = replication_origin_params.get_origin_interval()
    config.replication_terminus = (
            replication_origin_params.get_terminus_interval())
    config.ss_calculator_bin = settings.HYBRID_SS_MIN_BIN

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
