"""
Utility methods for creating JBrowse config files to allow data
to be viewed using JBrowse.
"""

import json
import os
import subprocess


from import_util import generate_fasta_from_genbank
from main.models import Dataset
from main.models import get_dataset_with_type
from main.models import ReferenceGenome
from settings import JBROWSE_BIN_PATH
from settings import JBROWSE_DATA_SYMLINK_PATH
from settings import JBROWSE_DATA_URL_ROOT
from settings import AWS_EXPECTED_BUCKET

# TODO: Figure out better place to put this.
# JBrowse requires the symlink path to exist. See settings.py
# comments for more info.

# assert os.path.exists(JBROWSE_DATA_SYMLINK_PATH), (
#         "%s does not exists. You may need to symlink it." %
#                 JBROWSE_DATA_SYMLINK_PATH)



def prepare_jbrowse_ref_sequence(reference_genome, **kwargs):
    """Prepare the reference sequence and place it in the ref_genome dir.

    This implicitly creates the config directory structure for this reference
    genome. Tracks added in the future are added relative to this reference
    genome.

    The implemenation of this method is a light wrapper around
    jbrowse/bin/prepare-refseqs.pl.
    """
    PREPARE_REFSEQS_BIN = os.path.join(JBROWSE_BIN_PATH, 'prepare-refseqs.pl')

    # First ensure that the reference genome exists. If it fails, try to
    # convert from genbank, then give up.
    reference_fasta = reference_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Next, ensure that the jbrowse directory exists.
    reference_genome.ensure_jbrowse_dir()
    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    # Now run prepare-refseqs.pl to get the ReferenceGenome in.
    subprocess.call([
        PREPARE_REFSEQS_BIN,
        '--fasta', reference_fasta,
        '--out', jbrowse_path,
    ])


def add_bam_file_track(reference_genome, sample_alignment, alignment_type):
    """Update the JBrowse track config file, trackList.json, for this
    ReferenceGenome with a track for the given sample_alignment and alignment_type.
    """
    ADD_TRACK_BIN = os.path.join(JBROWSE_BIN_PATH, 'add-track-json.pl')

    reference_genome = sample_alignment.alignment_group.reference_genome

    # Get the bam file location from the the Dataset of the genome
    # keyed by the alignment_type.
    bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)

    # Figure out the url that JBrowse would use to show the data, e.g.:
    #     /jbrowse/gd_data/projects/58a62c7d/genomes/8dc829ec/align.bam
    # urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
    #         bam_dataset.filesystem_location)

    # NOTE: We should construct bam file urls using project.get_client_jbrowse_link()
    #       rather than checking S3 flag here. 
    if reference_genome.project.is_s3_backed():
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % AWS_EXPECTED_BUCKET,
            bam_dataset.filesystem_location.strip("/jbrowse"))
    else:
        urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
            bam_dataset.filesystem_location)

    # Generic label for now.
    # TODO: Is there a better way to come up with a label?
    label = str(sample_alignment.experiment_sample.uid) + '_' + alignment_type

    # Build the JSON object.
    raw_dict_obj = {
        'storeClass': 'JBrowse/Store/SeqFeature/BAM',
        'urlTemplate': urlTemplate,
        'label': label,
        'type': 'JBrowse/View/Track/Alignments2',
        'key': label,
        'style' : {
            'className': 'alignment',
            'arrowheadClass': 'arrowhead',
            'labelScale': 100
        }
    }
    json_obj_string = json.dumps(raw_dict_obj)
    _update_ref_genome_track_list(reference_genome, json_obj_string)

    # Also add a snp coverage track.
    snp_coverage_label = label + '_coverage'
    coverage_raw_dict_obj = {
        'storeClass': 'JBrowse/Store/SeqFeature/BAM',
        'urlTemplate': urlTemplate,
        'label': snp_coverage_label,
        'type': 'JBrowse/View/Track/SNPCoverage',
        'key': snp_coverage_label,
    }
    coverage_json_obj_string = json.dumps(coverage_raw_dict_obj)
    _update_ref_genome_track_list(reference_genome, coverage_json_obj_string)


def _update_ref_genome_track_list(reference_genome, json_obj_string):
    """Helper method for updating the track list file for a ReferenceGenome.

    Args:
        reference_genome: A ReferenceGenome object.
        json_obj_string: The stringified json containing the track data.
    """
    ADD_TRACK_BIN = os.path.join(JBROWSE_BIN_PATH, 'add-track-json.pl')

    # Validate inputs.
    assert isinstance(reference_genome, ReferenceGenome)
    assert isinstance(json_obj_string, str)

    # The tracklist output file should be located relative to the
    # ReferenceGenome data dir.
    track_list_file = os.path.join(
            reference_genome.get_jbrowse_directory_path(), 'trackList.json')

    # Use JBrowse's utility script for adding a track.
    # The utility doesn't really do anything special beyond adding the json
    # string in the correct place in the file, but might as well use it since
    # it exists.
    echo_process = subprocess.Popen([
        'echo',
        json_obj_string
    ], stdout=subprocess.PIPE)

    subprocess.call([
        ADD_TRACK_BIN,
        track_list_file
    ], stdin=echo_process.stdout)
