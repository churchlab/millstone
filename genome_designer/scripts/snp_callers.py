"""
Functions for calling SNPs.
"""


def run_snp_calling_pipeline(alignment_group):
    """Calls SNPs for all of the alignments in the alignment_group.
    """
    # First make sure all of the samples have complete alignments.
    print 'Run SNP calling'

    # Validatate that all the necessary datasets in the alignment group
    # are ready.
    sample_alignment_list = (
            alignment_group.experimentsampletoalignment_set.all())
    for sample_alignment in sample_alignment_list:
        print sample_alignment
