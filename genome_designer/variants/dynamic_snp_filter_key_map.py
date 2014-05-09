"""
Script that updates the dynamic SNP Filter Key Map.

In our setup, we may want to query against data for an entire Variant, data
captured in VaiantCallerCommonData, or data describing the relationship between
a single Variant and ExperimentSample, captured in VariantEvidence.

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

import vcf

from main.exceptions import InputError
from main.models import ReferenceGenome
from variants.filter_key_map_constants import EXPERIMENT_SAMPLE_HARD_CODED
from variants.filter_key_map_constants import MAP_KEY__ALTERNATE
from variants.filter_key_map_constants import MAP_KEY__COMMON_DATA
from variants.filter_key_map_constants import MAP_KEY__EVIDENCE
from variants.filter_key_map_constants import MAP_KEY__EXPERIMENT_SAMPLE
from variants.filter_key_map_constants import SNP_CALLER_COMMON_DATA_HARD_CODED
from variants.filter_key_map_constants import SNP_EVIDENCE_HARD_CODED
from variants.filter_key_map_constants import SNP_VARIANT_HARD_CODED


def initialize_filter_key_map():
    """Initialize the filter key map with hard-coded fields. This should
    run on a signal of new ref genome creation.
    """

    variant_key_map = copy.deepcopy({
        MAP_KEY__ALTERNATE:
                SNP_VARIANT_HARD_CODED,
        MAP_KEY__COMMON_DATA:
                SNP_CALLER_COMMON_DATA_HARD_CODED,
        MAP_KEY__EVIDENCE:
                SNP_EVIDENCE_HARD_CODED,
        MAP_KEY__EXPERIMENT_SAMPLE:
                EXPERIMENT_SAMPLE_HARD_CODED
    })

    _assert_unique_keys(variant_key_map)

    return variant_key_map


def update_sample_filter_key_map(ref_genome, experiment_sample):
    """
    Updates a reference genome's variant key map dictionary with
    (potentially) new fields from sample metadata.

    Every time a new alignment happens, we need to update this
    with potentially new sample metadata keys.

    Every sample field gets uppercase_underscore()ed and
    prefixed with SAMPLE_ for filtering purposes.

    #TODO: deal with casting. for now, all are strings.
    """

    # Get an updated representation of the reference genome, in case
    # other tasks have already updated the filter_key_map
    # https://github.com/churchlab/millstone/issues/254
    ref_genome = ReferenceGenome.objects.get(uid=ref_genome.uid)

    if not any(ref_genome.variant_key_map):
        print ref_genome.variant_key_map
        raise ValueError(
                'Variant Key Map was not properly initialized for %s' % (
                        ref_genome.uid))

    sample_map = copy.deepcopy(
            ref_genome.variant_key_map.get(
            MAP_KEY__EXPERIMENT_SAMPLE, {}))

    for key, value in experiment_sample.data.iteritems():
        inner_map = {}
        inner_map['type'] = 'String'
        inner_map['num'] = 1

        sample_map[key]=inner_map

    ref_genome.variant_key_map[MAP_KEY__EXPERIMENT_SAMPLE] = sample_map

    _assert_unique_keys(ref_genome.variant_key_map)
    ref_genome.save(update_fields=['variant_key_map'])


def update_filter_key_map(ref_genome, source_vcf):
    """Updates a reference genome's variant key map dictionary with
    (potentially) new fields from a new VCF.

    The strategy is to use pyvcf to do the heavy-lifting of parsing the vcf
    header and get the types for each key out of there.

    We store a map of these types in a ReferenceGenome.variant_key_map,
    whose underlying representation is a Postgres json column.

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

    # Get an updated representation of the reference genome, in case
    # other tasks have already updated the filter_key_map
    # https://github.com/churchlab/millstone/issues/254
    ref_genome = ReferenceGenome.objects.get(uid=ref_genome.uid)

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

    common_data_map = ref_genome.variant_key_map.get(MAP_KEY__COMMON_DATA, {})
    alternate_map = ref_genome.variant_key_map.get(MAP_KEY__ALTERNATE, {})
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
    ref_genome.variant_key_map.get(MAP_KEY__COMMON_DATA, {}).update(
            common_data_map)
    ref_genome.variant_key_map.get(MAP_KEY__ALTERNATE, {}).update(
            alternate_map)

    # Update all of the per-sample fields.
    evidence_data_map = ref_genome.variant_key_map.get(MAP_KEY__EVIDENCE, {})
    for orig_key, value in vcf_reader.formats.iteritems():
        key = orig_key
        inner_map = {}
        inner_map['type'] = value.type
        inner_map['num'] = value.num
        evidence_data_map[key] = inner_map

    ref_genome.variant_key_map.get(MAP_KEY__EVIDENCE, {}).update(
            evidence_data_map)

    _assert_unique_keys(ref_genome.variant_key_map)

    ref_genome.save(update_fields=['variant_key_map'])


def _assert_unique_keys(variant_key_map):
    """Checks that the keys are unique across different submaps.

    Other parts of the code make the assumption that they subkeys are unique.
    I'm open to changing this requirement as long as we do it safely.
    """
    field_name_sets = [(submap, set(variant_key_map[submap].keys()))
            for submap in variant_key_map.keys()]
    # Check all pairs for intersections.
    for i in range(len(field_name_sets[1])):
        for j in range(len(field_name_sets[1])):
            submap_i = field_name_sets[i][0]
            submap_j = field_name_sets[j][0]

            if i <= j:
                continue
            assert not (field_name_sets[i][1] & field_name_sets[j][1]), (
                'Duplicate SNP filter keys between submaps {submap_i}'
                ' and {submap_j}. Keys in commmon: {common_keys}'.format(
                    submap_i= submap_i,
                    submap_j = submap_j,
                    common_keys = ' '.join(
                            field_name_sets[i][1] & field_name_sets[j][1])
                )
            )
