"""Wrapper for running Freebayes.
"""

import subprocess

from django.conf import settings

from main.model_utils import get_dataset_with_type


def run_freebayes(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, **kwargs):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.

    Returns:
        Boolean, True if successfully made it to the end, else False.
    """
    bam_files = [
            get_dataset_with_type(sa, alignment_type).get_absolute_location()
            for sa in sample_alignments]

    # Build up the bam part of the freebayes binary call.
    bam_part = []
    for bam_file in bam_files:
        bam_part.append('--bam')
        bam_part.append(bam_file)

    # Determine alignment ploidy (haploid or diploid).
    alignment_group = sample_alignments[0].alignment_group
    if alignment_group.alignment_options['call_as_haploid']:
        alignment_ploidy = 1
    else:
        alignment_ploidy = 2

    other_args_part = [
        '--fasta-reference', fasta_ref,
        '--pvar', '0.001',
        '--ploidy', str(alignment_ploidy),
        '--min-alternate-fraction', '.3',
        '--hwe-priors-off',
        # '--binomial-obs-priors-off',
        '--use-mapping-quality',
        '--min-base-quality', '25',
        '--min-mapping-quality', '30'
    ]

    # Build the full command and execute it for all bam files at once.
    full_command = (
            ['%s/freebayes/freebayes' % settings.TOOLS_DIR] +
            bam_part +
            other_args_part)

    with open(vcf_output_filename, 'w') as fh:
        subprocess.check_call(full_command, stdout=fh)

    return True # success
