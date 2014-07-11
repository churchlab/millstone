"""Script to generate coverage data.
"""

import os
import subprocess

from django.conf import settings

from main.models import get_dataset_with_type
from main.models import AlignmentGroup
from main.models import Dataset
from utils import generate_safe_filename_prefix_from_label


def analyze_coverage(sample_alignment, output_dir):
    ref_genome_fasta_location = get_dataset_with_type(
            sample_alignment.alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    input_bam_file = sample_alignment.dataset_set.get(
            type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()

    output_filename = generate_safe_filename_prefix_from_label(
        sample_alignment.experiment_sample.label + '_' +
        sample_alignment.uid) + '.coverage'
    output_path = os.path.join(output_dir, output_filename)

    with open(output_path, 'w') as fh:
        p_mpileup = subprocess.Popen([
                '%s/samtools/samtools' % settings.TOOLS_DIR,
                'mpileup',
                '-f', ref_genome_fasta_location,
                input_bam_file
        ], stdout=subprocess.PIPE)

        subprocess.check_call([
                'cut',
                '-f',
                '-4'
        ], stdin=p_mpileup.stdout, stdout=fh)
