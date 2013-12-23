"""
Methods for reading vcf files and importing it in our database representation.

We leverage pyvcf as much as possible.
"""

import pickle

from django.db import transaction
import vcf

from main.models import get_dataset_with_type
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantAlternate
from main.models import VariantEvidence
from scripts.dynamic_snp_filter_key_map import update_filter_key_map


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


@transaction.commit_on_success
def parse_vcf(vcf_dataset, reference_genome):
    """Parses the VCF and creates Variant models relative to ReferenceGenome.
    """
    query_cache = QueryCache()

    with open(vcf_dataset.get_absolute_location()) as fh:
        vcf_reader = vcf.Reader(fh)

        # First, update the reference_genome's key list with any new
        # keys from this VCF.
        update_filter_key_map(reference_genome, vcf_reader)

        for record_idx, record in enumerate(vcf_reader):
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

def get_or_create_variant(reference_genome, vcf_record, vcf_dataset,
        query_cache):
    """Create a variant if it doesn't already exist, along with
    its child VariantCallerCommonData object and VariantEvidencye
    and VariantAlternate children.

    Right now this assumes we are always using Freebayes for alignment.

    Args:
        reference_genome: The ReferenceGenome.
        raw_data_dict: Data dictionary with at least the following keys:
            * CHROM
            * POS
            * REF
            * ALT

        Also go through all per-alt keys and add them as a json field
        to the VariantAlternate object.

    Returns:
        Tuple (Variant, List<VariantAlt>)
    """

    # Build a dictionary of data for this record.
    raw_data_dict = extract_raw_data_dict(vcf_record)

    # Extract the relevant fields from the common data object.
    type = 'UNKNOWN' # TODO: Get this from vcf data?
    chromosome = pickle.loads(str(raw_data_dict.pop('CHROM')))
    position = int(pickle.loads(str(raw_data_dict.pop('POS'))))
    ref_value = pickle.loads(str(raw_data_dict.pop('REF')))
    alt_values = pickle.loads(str(raw_data_dict.pop('ALT')))

    # Try to find an existing one, or create it.
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

        # Unpickle the list and then repickle the single index needed
        alt_data = dict([
                (k, pickle.dumps(pickle.loads(raw_data_dict[k])[alt_idx])) 
                for k in raw_alt_keys])

        var_alt, var_created = VariantAlternate.objects.get_or_create(
                variant=variant,
                alt_value=alt_value)

        # If this is a new alternate, initialize the data dictionary
        if var_created: var_alt.data = {}
        # TODO: we are overwriting keys here; not sure this is desired behavior
        var_alt.data.update(alt_data)
        # Removing the ALT key since it is stored 
        var_alt.save()

        alts.append(var_alt)

    # Remove all per-alt keys from raw_data_dict before passing to VCC create
    [raw_data_dict.pop(k, None) for k in raw_alt_keys]

    # Create a common data object for this variant.
    common_data_obj, created = VariantCallerCommonData.objects.get_or_create(
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
        ve, created = VariantEvidence.objects.get_or_create(
                experiment_sample=sample_obj,
                variant_caller_common_data=common_data_obj,
                data=sample_data_dict)

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
