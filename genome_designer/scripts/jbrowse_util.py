"""
Utility methods for creating JBrowse config files to allow data
to be viewed using JBrowse.
"""

import os
import subprocess

from main.models import Dataset
from settings import JBROWSE_BIN_PATH
from settings import JBROwSE_DATA_SYMLINK_PATH


# JBrowse requires the symlink path to exist. See settings.py
# comments for more info.
assert os.path.exists(JBROwSE_DATA_SYMLINK_PATH), (
        "%s does not exists. You may need to symlink it." %
                JBROwSE_DATA_SYMLINK_PATH)


def prepare_reference_sequence(reference_genome):
    """Prepare the reference sequence and place it in the ref_genome dir.

    This implicitly creates the config directory structure for this reference
    genome. Tracks added in the future are added relative to this reference
    genome.

    The implemenation of this method is a light wrapper around
    jbrowse/bin/prepare-refseqs.pl.
    """
    PREPARE_REFSEQS_BIN = os.path.join(JBROWSE_BIN_PATH, 'prepare-refseqs.pl')

    # First ensure that the reference genome exists.
    reference_fasta = reference_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    assert reference_fasta is not None, "No reference soure fasta."

    # Next, ensure that the jbrowse directory exists.
    reference_genome.ensure_jbrowse_dir()
    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    # Now run prepare-refseqs.pl to get the ReferenceGenome in.
    subprocess.call([
        PREPARE_REFSEQS_BIN,
        '--fasta', reference_fasta,
        '--out', jbrowse_path,
    ])
