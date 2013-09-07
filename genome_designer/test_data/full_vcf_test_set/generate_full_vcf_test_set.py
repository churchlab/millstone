"""
Script for generating the test set.

This document describes how this test test was generated.

1) Select a region of the MG1655 genome to excise.
"""

import copy
import random

from Bio import SeqIO
import vcf

import simNGS_util


# Portion of MG1655 Genbank of size ~5.5 kB
EXCISED_GENBANK = 'mg1655_tolC_through_zupT.gb'

TEMPLATE_VCF = 'template.vcf'
VCF_TEMPLATE_READER = vcf.Reader(TEMPLATE_VCF)

SAMPLE_FASTA_ROOT = 'sample'

DESIGNED_SNP_VCF = 'designed_snps.vcf'



# If we do a SNP every 100 bases, that's 50 SNPs.
# We'll then do 20 designed SNPs and 20 SNPs per sample so we should get
# fairly interesting overlaps.
TOTAL_SNPS = 50

NUM_IN_CDS = 45

NUM_OTHER = TOTAL_SNPS - NUM_IN_CDS


# We'll create this many genomes.
NUM_SAMPLES = 6


def is_position_in_coding_feature(position, cds_features):
    """Checks whether the given position lies inside of a coding feature
    in the given genome record.
    """
    for feature in cds_features:
        if (feature.location.start <= position and
                position < feature.location.end):
            return True
    return False

BASE_OPTIONS = ['A', 'T', 'G', 'C']

def choose_alt(ref):
    """Returns a random base that is not ref.
    """
    alt = ref
    while alt == ref:
        alt = random.choice(BASE_OPTIONS)
    return alt


def get_subset_of_snps(all_snps, subset_size):
    all_snp_positions = all_snps.keys()
    subset = {}
    while len(subset) < subset_size:
        pos = random.choice(all_snp_positions)
        if pos in subset:
            continue
        subset[pos] = all_snps[pos]
    return subset


def create_vcf_for_subset(subset, out_path):
    with open(out_path, 'w') as designed_fh:
        writer = vcf.Writer(designed_fh, VCF_TEMPLATE_READER,
                lineterminator='\n')
        for pos, value_dict in subset.iteritems():
            writer.write_record(vcf.model._Record(
                    'Chromosome', # CHROM
                    pos, # POS
                    None, # ID
                    value_dict['ref'], # REF
                    value_dict['alt'], # ALT
                    None, # QUAL
                    None, # FILTER
                    None, # INFO
                    None, # FORMAT
                    None, # sample_indexes
                    samples=None))


def main():
    seq_record = SeqIO.read(EXCISED_GENBANK, 'genbank')
    cds_features = [f for f in seq_record.features if f.type == 'CDS']

    # Generate all possible SNPs to sample from. Store them in a dictionary
    # keyed by position so we can easily deal with lookups and avoiding
    # duplicates as needed below.
    all_snps = {}

    len_seq_record = len(seq_record)

    # Select random positions for SNPs, respecting the distribution
    # set above by the NUM_IN_CDS vs TOTAL_SNPS constants.
    # NOTE: These SNP positions are pythonic. We have to update them when
    # writing them out in vcf format below.
    num_in_cds = 0
    num_other = 0
    while num_in_cds < NUM_IN_CDS or num_other < NUM_OTHER:
        position = random.randint(0, len_seq_record - 1)
        if position in all_snps:
            continue

        in_cds_feature = is_position_in_coding_feature(position, cds_features)
        do_add_position = False
        if in_cds_feature and num_in_cds < NUM_IN_CDS:
            do_add_position = True
            num_in_cds += 1
        elif not in_cds_feature and num_other < NUM_OTHER:
            do_add_position = True
            num_other += 1

        if do_add_position:
            ref = seq_record.seq[position]
            alt = choose_alt(ref)
            all_snps[position] = {
                'ref': ref,
                'alt': [alt]
            }

    assert len(all_snps) == TOTAL_SNPS, "Didn't get all the SNPs we expected."

    # Now select a subset of these SNPS to serve as designed.
    designed_snps = get_subset_of_snps(all_snps, 20)
    create_vcf_for_subset(designed_snps, DESIGNED_SNP_VCF)

    # Now create the samples.
    for sample_num in range(NUM_SAMPLES):
        sample_name = SAMPLE_FASTA_ROOT + str(sample_num)
        sample_record = copy.deepcopy(seq_record)
        sample_record.id = sample_name

        # Grab a subset of SNPs.
        sample_snps = get_subset_of_snps(all_snps, 20)

        # Introduce the mutations.
        for position, value_dict in sample_snps.iteritems():
            sample_record.seq = (
                    sample_record.seq[:position] +
                    value_dict['alt'][0] + 
                    sample_record.seq[position + 1:])
        assert len(sample_record) == len(seq_record), (
                "For now we are only doing mutations.")

        # Write out the sample fasta.
        sample_output = sample_name + '.fa'
        with open(sample_output, 'w') as out_fh:
            SeqIO.write(sample_record, out_fh, 'fasta')

        # Generate fake reads using simNGS.
        simLibrary_fasta = sample_name + '.simLibrary.fa'
        print sample_output, simLibrary_fasta
        simNGS_util.run_simLibrary(sample_output, simLibrary_fasta)

        # Generate reads using simNGS.
        output_fq = sample_name + '.simLibrary.fq'
        simNGS_util.run_paired_simNGS(simLibrary_fasta, output_fq)


if __name__ == '__main__':
    main()
