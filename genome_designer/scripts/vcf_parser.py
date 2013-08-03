"""
Methods for reading vcf files and importing it in our database representation.

We leverage pyvcf as much as possible.
"""

import pickle

import vcf

from main.models import get_dataset_with_type
from main.models import Variant
from main.models import VariantCallerCommonData


def parse_alignment_group_vcf(alignment_group, vcf_dataset_type):
    """Parses the VCF associated with the AlignmentGroup and saves data there.
    """
    reference_genome = alignment_group.reference_genome

    vcf_dataset = get_dataset_with_type(alignment_group, vcf_dataset_type)

    with open(vcf_dataset.get_absolute_location()) as fh:
        vcf_reader = vcf.Reader(fh)
        for record in vcf_reader:
            # Create a common data object for this variant.
            common_data_obj = VariantCallerCommonData.objects.create(
                    reference_genome=reference_genome,
                    source_dataset=vcf_dataset,
            )
            # Extracting the data is somewhat of a process so we do
            # so in this function.
            _populate_common_data_from_vcf_record(common_data_obj, record)

    aggregate_variants(alignment_group)


def _populate_common_data_from_vcf_record(common_data_obj, vcf_record):
    """Get the desired data out of the vcf format and write it to
    common_data_obj.
    """
    # Keys that we need to do extra work with in order to copy.
    MANUAL_KEYS = ['INFO', 'samples']

    data_dict = {}

    # Copy over non-manual keys
    for key, value in vcf_record.__dict__.iteritems():
        if not key in MANUAL_KEYS:
            # For now, we pickle (and later un-pickle values) because
            # they come in different types, and we haven't figured out a good
            # way to deal with this quite yet.
            data_dict[str(key)] = pickle.dumps(value)
    common_data_obj.data = data_dict

    # Populate 'INFO'
    _populate_common_data_info(common_data_obj, vcf_record)

    # Commit the transaction.
    common_data_obj.save()


def _populate_common_data_info(common_data_obj, vcf_record):
    """Parses the vcf record INFO field and saves the result into the
    common data object data field.
    """
    data_dict = common_data_obj.data
    for key, value in vcf_record.INFO.iteritems():
        effective_key = 'INFO_' + key
        data_dict[effective_key] = pickle.dumps(value)
    common_data_obj.data = data_dict


def aggregate_variants(alignment_group):
    """Determine the unique set of Variants for this alignment group.
    """
    common_data_obj_list = VariantCallerCommonData.objects.filter(
            reference_genome=alignment_group.reference_genome)
    for common_data_obj in common_data_obj_list:
        get_or_create_variant(alignment_group, common_data_obj)


def get_or_create_variant(alignment_group, common_data_obj):
    """Create a variant if it doesn't already exist.

    Right now this assumes we are always using Freebayes for alignment.
    """
    raw_data_dict = common_data_obj.data

    # Extract the relevant fields from the common data object.
    type = 'UNKNOWN' # TODO: Get this from vcf data?
    chromosome = pickle.loads(str(raw_data_dict['CHROM']))
    position = int(pickle.loads(str(raw_data_dict['POS'])))
    ref_value = pickle.loads(str(raw_data_dict['REF']))
    alt_value = pickle.loads(str(raw_data_dict['ALT']))


    # Try to find an existing one, or create it.
    variant, created = Variant.objects.get_or_create(
            reference_genome=alignment_group.reference_genome,
            type=type,
            chromosome=chromosome,
            position=position,
            ref_value=ref_value,
            alt_value=alt_value,
    )

    # Link to this common_data object.
    common_data_obj.variant = variant
    common_data_obj.save()

    return variant
