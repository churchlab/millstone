"""
Script that updates the dynamic SNP Filter Key Map.

In our setup, we may want to query against data for an entire Variant, data
captured in VaiantCallerCommonData, or data describing the relationship between
a single Variant and ExperimentSample, captured in VariantEvidence.


Note that the fields are stored in a key-value
field in the database, and because the values may be of different types,
we just pickle the pyvcf object to store it.

Previously, we were using a static map that we defined once, in
generate_filter_key_map.py. Now, because every VCF we import might have
different fields, and because SnpEFF adds fields, we need to dynamically
update all possible VCF fields per ReferenceGenome object. These fields
will be stored in a JSONField in the ReferenceGenome called
variant_key_map.

The maps defined here specify:
    * the valid keys that can be filtered agains
    * which data object they are located in, SNPCallerCommonData or SNPEvidence
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
    'ALT': {'type': 'String', 'num': -2}
}


SNP_EVIDENCE_HARD_CODED = {
    'gt_type': {'type': 'Integer', 'num': 1},
    'is_het': {'type': 'Boolean', 'num': 1}
}


def initialize_filter_key_map(ref_genome):
    """Initialize the filter key map with hard-coded fields. This should
    run on a signal of new ref genome creation.
    """
    ref_genome.variant_key_map = {
        'snp_caller_common_data': copy.deepcopy(
                SNP_CALLER_COMMON_DATA_HARD_CODED),
        'snp_evidence_data': copy.deepcopy(
                SNP_EVIDENCE_HARD_CODED)
    }
    ref_genome.save()


def update_filter_key_map(ref_genome, source_vcf):
    """Updates a reference genome's variant key map dictionary with
    (potentially) new fields from a new VCF.

    The strategy is to use pyvcf to do the heavy-lifting of parsing the vcf
    header and get the types for each key out of there.

    We store a map of these types in a JSONField object per reference-genome.

    TODO: Handle potential INFO ID name collisions. Probably won't happen. ;)
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

    # Write the map with the keys for SNPCallerCommonData
    snp_caller_common_data_map = ref_genome.variant_key_map[
            'snp_caller_common_data']
    for orig_key, value in vcf_reader.infos.iteritems():
        key = 'INFO_' + orig_key
        inner_map = {}
        inner_map['type'] = value.type
        inner_map['num'] = value.num
        snp_caller_common_data_map[key] = inner_map

    # Write the map with the keys for SNPEvidence
    snp_evidence_data_map = ref_genome.variant_key_map[
            'snp_evidence_data']
    for orig_key, value in vcf_reader.formats.iteritems():
        key = orig_key
        inner_map = {}
        inner_map['type'] = value.type
        inner_map['num'] = value.num
        snp_evidence_data_map[key] = inner_map

    # Test the generated file by running imports.
    ref_genome.variant_key_map['snp_caller_common_data'].update(
            snp_caller_common_data_map)
    ref_genome.variant_key_map['snp_evidence_data'].update(
           snp_evidence_data_map)
    ref_genome.save()
