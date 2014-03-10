"""
Methods for reading vcf files and importing it in our database representation.

We leverage pyvcf as much as possible.
"""

from django.db import transaction
import vcf

from main.models import get_dataset_with_type
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantAlternate
from main.models import VariantEvidence
from scripts.dynamic_snp_filter_key_map import update_filter_key_map


SV_TYPES = {
    'DEL': 'DELETION',
    'DUP': 'DUPLICATION',
    'INV': 'INVERSION'
}

UNKNOWN_VARIANT_TYPE = 'unknown'

IGNORE_VCF_RECORD_KEYS = [
    # Bug in pyvcf where some alleles are some sort of object rather than str.
    'alleles'
]


class QueryCache(object):
    """Manual caching for queries to avoid excessive db calls.

    NOTE: Only trying this with one model so far. Haven't determined how
    useful this really is.
    """
    def __init__(self):
        self.uid_to_experiment_sample_map = {}


def parse_alignment_group_vcf(alignment_group, vcf_dataset_type):
    """Parses the VCF associated with the AlignmentGroup and saves data there.
    """
    vcf_dataset = get_dataset_with_type(alignment_group, vcf_dataset_type)
    parse_vcf(vcf_dataset, alignment_group.reference_genome)


# @transaction.commit_on_success
def parse_vcf(vcf_dataset, reference_genome):
    """Parses the VCF and creates Variant models relative to ReferenceGenome.
    """
    # This helper object will help prevent repeated calls to the database.
    # We'll use it at least for ExperimentSamples.
    query_cache = QueryCache()

    # First count the number of records to give helpful status debug output.
    record_count = 0
    with open(vcf_dataset.get_absolute_location()) as fh:
        vcf_reader = vcf.Reader(fh)
        for record in vcf_reader:
            record_count += 1

    # Now iterate through the vcf file again and parse the data.
    with open(vcf_dataset.get_absolute_location()) as fh:
        vcf_reader = vcf.Reader(fh)

        # First, update the reference_genome's key list with any new
        # keys from this VCF.
        update_filter_key_map(reference_genome, vcf_reader)

        for record_idx, record in enumerate(vcf_reader):
            print 'vcf_parser: Parsing %d out of %d' % (
                    record_idx + 1, record_count)
            # Make sure the QueryCache object has experiment samples populated.
            # Assumes every row has same samples. (Pretty sure this is true
            # for well-formatted vcf file.)
            if (len(query_cache.uid_to_experiment_sample_map) == 0 and
                    len(record.samples) > 0):
                for sample in record.samples:
                    sample_uid = sample.sample
                    query_cache.uid_to_experiment_sample_map[sample_uid] = (
                            ExperimentSample.objects.get(uid=sample_uid))

            # Get or create the Variant for this record. This step
            # also generates the alternate objects and assigns their
            # data fields as well.
            variant, alts = get_or_create_variant(reference_genome,
                    record, vcf_dataset, query_cache)


def extract_raw_data_dict(vcf_record):
    """Extract a dictionary of raw data from the Record.

    Returns:
        Dictionary representation of the record.
    """
    # Keys that we need to do extra work with in order to copy.
    MANUAL_KEYS = ['INFO', 'samples']

    data_dict = {}

    # Extract keys that do not need to be handled manually.
    for key, value in vcf_record.__dict__.iteritems():
        if key in IGNORE_VCF_RECORD_KEYS or key in MANUAL_KEYS:
            continue
        data_dict[str(key)] = value

    # The TYPE is just a property of the record object.
    # TODO: Do we care about the var_subtype()? (ts/tv/complex/sv type/etc?)
    # sv callers store the type in vcf_record.INFO['SVTYPE']
    if hasattr(vcf_record, 'INFO') and 'SVTYPE' in dict(vcf_record.INFO):
        data_dict['TYPE'] = SV_TYPES.get(dict(vcf_record.INFO)['SVTYPE'],
                UNKNOWN_VARIANT_TYPE)
    elif hasattr(vcf_record, 'var_type'):
        data_dict['TYPE'] = str(vcf_record.var_type)
    else:
        data_dict['TYPE'] = UNKNOWN_VARIANT_TYPE

    # Populate 'INFO'
    if hasattr(vcf_record, 'INFO'):
        populate_common_data_info(data_dict, vcf_record)

    return data_dict


def populate_common_data_info(data_dict, vcf_record):
    """Parses the vcf record INFO field and updates the data dict.
    """
    for key, value in vcf_record.INFO.iteritems():
        effective_key = 'INFO_' + key
        data_dict[effective_key] = value


def get_or_create_variant(reference_genome, vcf_record, vcf_dataset,
        query_cache=None):
    """Create a variant and its relations.

    A new Variant is only created if it doesn't exist already. The following
    relations are created:
        * VariantCallerCommonData
        * VariantEvidence
        * VariantAlternate

    Also go through all per-alt keys and add them as a json field
    to the VariantAlternate object.

    Right now this assumes we are always using Freebayes for alignment.

    Args:
        reference_genome: The ReferenceGenome.
        vcf_record: pyvcf Record object.
        vcf_dataset: Source Dataset for this data.
        query_cache: QueryCache helper object for making queries.

    Returns:
        Tuple (Variant, List<VariantAlt>)
    """
    # Build a dictionary of data for this record.
    raw_data_dict = extract_raw_data_dict(vcf_record)

    # Extract the REQUIRED fields from the common data object.
    type = str(raw_data_dict.pop('TYPE'))
    chromosome = raw_data_dict.pop('CHROM')
    position = int(raw_data_dict.pop('POS'))
    ref_value = raw_data_dict.pop('REF')
    alt_values = raw_data_dict.pop('ALT')

    # Try to find an existing Variant, or create it.
    variant, created = Variant.objects.get_or_create(
            reference_genome=reference_genome,
            type=type,
            chromosome=chromosome,
            position=position,
            ref_value=ref_value
    )

    alts = []
    all_alt_keys = reference_genome.get_variant_alternate_map().keys()
    raw_alt_keys = [k for k in raw_data_dict.keys() if k in all_alt_keys]

    for alt_idx, alt_value in enumerate(alt_values):

        # Grab the alt data for this alt index.
        alt_data = dict([(k, raw_data_dict[k][alt_idx]) for k in raw_alt_keys])

        var_alt, var_created = VariantAlternate.objects.get_or_create(
                variant=variant,
                alt_value=alt_value)

        # If this is a new alternate, initialize the data dictionary
        if var_created:
            var_alt.data = {}

        # TODO: We are overwriting keys here. Is this desired?
        var_alt.data.update(alt_data)
        var_alt.save()

        alts.append(var_alt)

    # Remove all per-alt keys from raw_data_dict before passing to VCC create.
    [raw_data_dict.pop(k, None) for k in raw_alt_keys]

    # Create a common data object for this variant.
    # NOTE: raw_data_dict only contains the values that were not popped until
    # this point.
    common_data_obj = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset=vcf_dataset,
            data=raw_data_dict
    )

    # TODO: What about num -2 objects? I'm really not excited about
    # creating a VariantGenotype object, nor do I think it will
    # be necessary, so skipping it for now, and keeping that data in
    # the VCC object.

    # Create a VariantEvidence object for each ExperimentSample.
    # NOTE: VariantEvidence are automatically linked to the correct
    #     VariantAlternate after they are created in
    #     main.signals.post_variant_evidence_create()
    for sample in vcf_record.samples:
        sample_uid = sample.sample
        sample_data_dict = extract_sample_data_dict(sample)
        if query_cache is not None:
            sample_obj = query_cache.uid_to_experiment_sample_map[sample_uid]
        else:
            sample_obj = ExperimentSample.objects.get(uid=sample_uid)
        VariantEvidence.objects.create(
                experiment_sample=sample_obj,
                variant_caller_common_data=common_data_obj,
                data=sample_data_dict)

    return (variant, alts)


def extract_sample_data_dict(s):
    """Extract sample data from the pyvcf _Call object.

    Args:
        s: pyvcf _Call object.

    Returns:
        A dictionary representing the object.
    """
    def _add_property(result_dict, s, key, eval_string):
        """Helper method to add keys. PyVCF is really buggy so we need to
        be extra paranoid when parsing it.
        """
        try:
            result_dict[key] = getattr(s, eval_string)
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
    )
    for key, eval_string in key_eval_string_pairs:
        _add_property(result, s, key, eval_string)

    # Add data fields in slightly different way.
    SAMPLE_DATA_FIELDS = [
        'AO',
        'DP',
        'GL',
        'GT',
        'QA',
        'QR',
        'RO',
    ]
    if hasattr(s, 'data'):
        for key in SAMPLE_DATA_FIELDS:
            if hasattr(s.data, key):
                result[key] = getattr(s.data, key)

    # TODO: Add support for SnpEff data.

    return result
