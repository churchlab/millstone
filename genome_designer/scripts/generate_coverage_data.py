"""Script to generate coverage data.
"""

import os
import subprocess

from main.models import get_dataset_with_type
from main.models import AlignmentGroup
from main.models import Dataset
from utils import generate_safe_filename_prefix_from_label
from utils.samtools_utils import run_mpileup


def analyze_coverage(sample_alignment, output_dir, coverage_only=True):
    input_ref_genome_fasta_location = get_dataset_with_type(
            sample_alignment.alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    input_bam_file = sample_alignment.dataset_set.get(
            type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()

    output_filename = generate_safe_filename_prefix_from_label(
        sample_alignment.experiment_sample.label + '_' +
        sample_alignment.uid) + '.coverage'
    output_path = os.path.join(output_dir, output_filename)

    #input_ref_genome_fasta_path = sample_alignment.alignment_group.reference_genome.dataset_set.get(type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    run_mpileup(input_bam_file, input_ref_genome_fasta_path, output_path,
            coverage_only=coverage_only)
