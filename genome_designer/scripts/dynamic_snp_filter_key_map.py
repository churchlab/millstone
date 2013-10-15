"""
Script that updates the dynamic SNP Filter Key Map.

In our setup, we may want to query against data for an entire Variant, data
captured in VaiantCallerCommonData, or data describing the relationship between
a single Variant and ExperimentSample, captured in VariantEvidence.


Note that the fields are stored in a key-value field in the database, and
because the values may be of different types, we just pickle the pyvcf object
to store it.

Previously, we were using a static map that we defined once, in
generate_filter_key_map.py. Now, because every VCF we import might have
different fields, and because SnpEFF adds fields, we need to dynamically
update all possible VCF fields per ReferenceGenome object. These fields will
be stored in a JSONField in the ReferenceGenome called variant_key_map.

The maps defined here specify:
    * the valid keys that can be filtered agains
    * which data object they are located in, SNPCallerCommonData, 
        SNP_Alternate_Data,
        and SNPEvidence
    * the type and number of the field (i.e. array or single, integer, float,
        or string)
    * the valid operations that can be performed on the types.
"""

import copy
import pprint

import vcf

SNP_CALLER_COMMON_DATA_HARD_CODED = {
    'CHROM': {'type': 'String', 'num': 1},
    'POS': {'type': 'Integer', 'num': 1},
    'REF': {'type': 'String', 'num': 1},
}

SNP_VARIANT_HARD_CODED = {
    'ALT': {'type': 'String', 'num': -1}
}

SNP_EVIDENCE_HARD_CODED = {
    'gt_type': {'type': 'Integer', 'num': 1},
    'is_het': {'type': 'Boolean', 'num': 1}
}

MAP_KEY__VARIANT = 'variant_data'

MAP_KEY__ALTERNATE = 'snp_alternate_data'

MAP_KEY__COMMON_DATA = 'snp_caller_common_data'

MAP_KEY__EVIDENCE = 'snp_evidence_data'


def initialize_filter_key_map(ref_genome):
    """Initialize the filter key map with hard-coded fields. This should
    run on a signal of new ref genome creation.
    """
    ref_genome.variant_key_map = {
        MAP_KEY__ALTERNATE: copy.deepcopy(
                SNP_VARIANT_HARD_CODED),
        MAP_KEY__COMMON_DATA: copy.deepcopy(
                SNP_CALLER_COMMON_DATA_HARD_CODED),
        MAP_KEY__EVIDENCE: copy.deepcopy(
                SNP_EVIDENCE_HARD_CODED)
    }
    ref_genome.save()


def update_filter_key_map(ref_genome, source_vcf):
    """Updates a reference genome's variant key map dictionary with
    (potentially) new fields from a new VCF.

    The strategy is to use pyvcf to do the heavy-lifting of parsing the vcf
    header and get the types for each key out of there.

    We store a map of these types in a JSONField object per reference-genome.

    TODO: Handle potential INFO ID name collisions (i.e. two different
        vcf files that have two non-identical ID fields with the same name.
        Probably won't happen. ;)

    Note about 'num': it's handled by pyvcf in a dict called `field_counts`, and
    'A' is stored as -1, which means that a single value is held for every
    ALTERATE allele. 'G' is stored as -2, which means there is a different value
    for every possible  genotype combination of alleles (which would be a choose
    n where n is the called  ploidy and a is the number of alleles). 

    All of the -1 fields are stored in a subdict corresponding to
    MAP_KEY__ALTERNATE, and the JSONField data in the VariantEvidence object
    stores the per-alt data after doing the equivalent of 'zip()ing' it per
    object.

    """

    #First try the source_vcf as a vcf file
    try:
        vcf_source_fh = open(source_vcf)
        vcf_reader = vcf.Reader(vcf_source_fh)
    except:
        #Then try the source_vcf as a vcf_reader
        try:
            assert hasattr(source_vcf, 'infos')
            assert hasattr(source_vcf, 'formats')
            vcf_reader = source_vcf
        except:
            raise InputError('Bad source_vcf arg: not filename or vcf_reader')

    common_data_map = ref_genome.variant_key_map[MAP_KEY__COMMON_DATA]
    alternate_map = ref_genome.variant_key_map[MAP_KEY__ALTERNATE]
    for orig_key, value in vcf_reader.infos.iteritems():
        key = 'INFO_' + orig_key
        inner_map = {}
        inner_map['type'] = value.type
        inner_map['num'] = value.num

        # If a field is per-alternate, then put it in a separate dictionary.
        if value.num == -1:
            alternate_map[key] = inner_map
        else:
            common_data_map[key] = inner_map

    # Update the reference genome's variant key maps with these (ostensibly)
    # new fields, overwriting previous data
    ref_genome.variant_key_map[MAP_KEY__COMMON_DATA].update(common_data_map)
    ref_genome.variant_key_map[MAP_KEY__ALTERNATE].update(alternate_map)

    # Update all of the per-sample fields.
    evidence_data_map = ref_genome.variant_key_map[MAP_KEY__EVIDENCE]
    for orig_key, value in vcf_reader.formats.iteritems():
        key = orig_key
        inner_map = {}
        inner_map['type'] = value.type
        inner_map['num'] = value.num
        evidence_data_map[key] = inner_map

    ref_genome.variant_key_map[MAP_KEY__EVIDENCE].update(evidence_data_map)

    ref_genome.save()
