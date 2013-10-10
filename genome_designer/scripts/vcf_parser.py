"""
Methods for reading vcf files and importing it in our database representation.

We leverage pyvcf as much as possible.
"""

import pickle
import vcf

from main.models import get_dataset_with_type
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantAlternate
from main.models import VariantEvidence
from scripts.dynamic_snp_filter_key_map import update_filter_key_map


def parse_alignment_group_vcf(alignment_group, vcf_dataset_type):
    """Parses the VCF associated with the AlignmentGroup and saves data there.
    """
    reference_genome = alignment_group.reference_genome

    vcf_dataset = get_dataset_with_type(alignment_group, vcf_dataset_type)

    with open(vcf_dataset.get_absolute_location()) as fh:
        vcf_reader = vcf.Reader(fh)
        for record in vcf_reader:

            # Build a dictionary of data for this record.
            raw_data_dict = extract_raw_data_dict(record)

            # Get or create the Variant for this record.
            variant, alts = get_or_create_variant(reference_genome, 
                    raw_data_dict)

            # Create a common data object for this variant.
            common_data_obj = VariantCallerCommonData.objects.create(
                    variant=variant,
                    source_dataset=vcf_dataset,
                    data=raw_data_dict
            )

            # Create a VariantEvidence object for each ExperimentSample.
            for sample in record.samples:
                sample_uid = sample.sample
                sample_data_dict = extract_sample_data_dict(sample)
                sample_obj = ExperimentSample.objects.get(uid=sample_uid)
                VariantEvidence.objects.create(
                        experiment_sample=sample_obj,
                        variant_caller_common_data=common_data_obj,
                        data=sample_data_dict)

        # Finally, update the reference_genome's key list with any new
        # keys from this VCF.
        update_filter_key_map(reference_genome, vcf_reader)


def extract_raw_data_dict(vcf_record):
    """Extract a dictionary of raw data from the Record.

    Returns:
        Dictionary representation of the record.
    """
    # Keys that we need to do extra work with in order to copy.
    MANUAL_KEYS = ['INFO', 'samples', 'ALT']

    data_dict = {}

    # Copy over non-manual keys
    for key, value in vcf_record.__dict__.iteritems():
        if not key in MANUAL_KEYS:
            # For now, we pickle (and later un-pickle values) because
            # they come in different types, and we haven't figured out a good
            # way to deal with this quite yet.
            data_dict[str(key)] = pickle.dumps(value)

    # The ALT key is a fancy _AltRecord object (or possibly subclassed 
    # as a substitution or an SV, so for now, just pass its __str__()
    data_dict['ALT'] = pickle.dumps(vcf_record.ALT)

    # The TYPE is just a property of the record object.
    # TODO: Do we care about the var_subtype()? (ts/tv/complex/sv type/etc?)
    if hasattr(vcf_record, 'var_type'):
        data_dict['TYPE'] = pickle.dumps(str(vcf_record.var_type))
    else:
        data_dict['TYPE'] =  'unknown'

    # Populate 'INFO'
    if hasattr(vcf_record, 'INFO'):
        populate_common_data_info(data_dict, vcf_record)

    return data_dict


def populate_common_data_info(data_dict, vcf_record):
    """Parses the vcf record INFO field and updates the data dict.
    """
    for key, value in vcf_record.INFO.iteritems():
        effective_key = 'INFO_' + key
        data_dict[effective_key] = pickle.dumps(value)

def get_or_create_variant(reference_genome, raw_data_dict):
    """Create a variant if it doesn't already exist.

    Right now this assumes we are always using Freebayes for alignment.

    Args:
        reference_genome: The ReferenceGenome.
        raw_data_dict: Data dictionary with at least the following keys:
            * CHROM
            * POS
            * REF
            * ALT
    """
    # Extract the relevant fields from the common data object.
    type = 'UNKNOWN' # TODO: Get this from vcf data?
    chromosome = pickle.loads(str(raw_data_dict['CHROM']))
    position = int(pickle.loads(str(raw_data_dict['POS'])))
    ref_value = pickle.loads(str(raw_data_dict['REF']))
    alt_values = pickle.loads(raw_data_dict['ALT'])

    # Try to find an existing one, or create it.
    variant, created = Variant.objects.get_or_create(
            reference_genome=reference_genome,
            type=type,
            chromosome=chromosome,
            position=position,
            ref_value=ref_value
    )

    alts = []

    for alt_value in alt_values:
        var_alt, va_created = VariantAlternate.objects.get_or_create(
            variant=variant,
            alt_value=alt_value)
        alts.append(var_alt)

    return (variant, alts)

def extract_sample_data_dict(s):
    """Manually serializes a pyvcf _Call object because their internal use of
    __slots__ breaks python pickle.

    Args:
        pyvcf _Call object.

    Returns:
        A dictionary representing the object.
    """
    def _add_property(result_dict, s, key, eval_string):
        """Helper method to add keys. PyVCF is really buggy so we need to
        be extra paranoid when parsing it.
        """
        try:
            result_dict[key] = eval('pickle.dumps(s.' + eval_string + ')')
        except AttributeError:
            result_dict[key] = None

    result = {}

    # The main properties we'll want to query across.
    key_eval_string_pairs = (
        ('sample', 'sample'),
        ('called', 'called'),
        ('gt_bases', 'gt_bases'),
        ('gt_nums', 'gt_nums'),
        ('gt_type', 'gt_type'),
        ('is_het', 'is_het'),
        ('is_variant', 'is_variant'),
        ('phased', 'phased'),
        ('AO', 'data.AO'),
        ('DP', 'data.DP'),
        ('GL', 'data.GL'),
        ('GT', 'data.GT'),
        ('QA', 'data.QA'),
        ('QR', 'data.QR'),
        ('RO', 'data.RO')
    )

    # TODO: Add support for SnpEff data. I don't remember why I was hard-coding
    # the above fields before, but it would be nice to avoid hard-coding if
    # possible.

    for key, eval_string in key_eval_string_pairs:
        _add_property(result, s, key, eval_string)
    return result
