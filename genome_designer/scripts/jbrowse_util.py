"""
Utility methods for creating JBrowse config files to allow data
to be viewed using JBrowse.
"""

import json
import os
import subprocess

from main.models import Dataset
from main.models import get_dataset_with_type
from main.models import ReferenceGenome
from settings import JBROWSE_BIN_PATH
from settings import JBROWSE_DATA_URL_ROOT
from settings import S3_BUCKET
from settings import TOOLS_DIR

TABIX_BINARY = '%s/tabix/tabix' % TOOLS_DIR

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
    reference_fasta = get_dataset_with_type(
            reference_genome,
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

def add_genbank_file_track(reference_genome, **kwargs):
    """
    Jbrowse has the ability to make tracks out of genbank files. This
    takes the genbank file from a reference_genome object and attempts to 
    make such a track and then add it to the track list.
    """
    FLATFILE_TRACK_BIN = os.path.join(JBROWSE_BIN_PATH, 'flatfile-to-json.pl')

    reference_gbk = get_dataset_with_type(
            reference_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    genbank_json_command = [
        FLATFILE_TRACK_BIN,
        '--gbk', reference_gbk,
        '--out', jbrowse_path,
        #'--type', 'CDS',
        '--autocomplete','all',
        '--trackLabel','gbk',
        '--key',"GenBank CDS",
        '--trackType',"JBrowse/View/Track/HTMLFeatures",
        #'--getSubfeatures',
        #'--className','transcript',
        '--subfeatureClasses', "{\"CDS\":\"transcript-CDS\"}"
    ]

    print ' '.join(genbank_json_command)

    subprocess.call(genbank_json_command)

def add_vcf_track(reference_genome, alignment_group, vcf_dataset_type):
    """
    From

    JBrowse Docs:
        http://gmod.org/wiki/JBrowse_Configuration_Guide
                #Example_VCF-based_Variant_Track_Configuration
    """
    # Get the vcf file location from the the Dataset of the genome
    # keyed by the alignment_type.
    vcf_dataset = get_dataset_with_type(alignment_group, 
            vcf_dataset_type)

    vcf_dataset = _vcf_to_vcftabix(vcf_dataset)

    if reference_genome.project.is_s3_backed():
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            vcf_dataset.filesystem_location.strip("/jbrowse"))
        urlTemplate_idx = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            vcf_dataset.filesystem_idx_location.strip("/jbrowse"))
    else:
        urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
            vcf_dataset.filesystem_location)
        urlTemplate_idx = os.path.join(JBROWSE_DATA_URL_ROOT,
            vcf_dataset.filesystem_idx_location)

    # TODO: This jbrowse vcf label really need to be more human readable.
    label = str(alignment_group.uid) + '_' + vcf_dataset.type

    # Build the JSON object.
    raw_dict_obj = {
        "label"         : label,
        "key"           : "%s SNVs" % alignment_group.label,
        "storeClass"    : "JBrowse/Store/SeqFeature/VCFTabix",
        "urlTemplate"   : urlTemplate,
        "tbiUrlTemplate": urlTemplate_idx,
        "type"          : "JBrowse/View/Track/HTMLVariants"
    }

    json_obj_string = json.dumps(raw_dict_obj)
    _update_ref_genome_track_list(reference_genome, json_obj_string)

def _vcf_to_vcftabix(vcf_dataset):
    """
    Generate a compressed version of a vcf file and index it with 
    samtools tabix. Add a new dataset model instance for this compressed
    version, with the same related objects. Flag is as compressed,
    indexed, etc. 
    """

    # 1. Check if dataset is compressed. If not, then grab the compressed
    # version or make a compressed version.
    if not vcf_dataset.is_compressed():
        # Check for existing compressed version using related model.
        # Assume that the first model will do. 
        related_model = vcf_dataset.get_related_model_set().all()[0]
        compressed_dataset = get_dataset_with_type(
                entity= related_model, 
                type= vcf_dataset.type, 
                compressed=True)
        # If there is no compressed dataset, then make it
        if compressed_dataset is None:
            compressed_dataset = vcf_dataset.make_compressed('.bgz')
    else:
        compressed_dataset = vcf_dataset

    if compressed_dataset.filesystem_idx_location == '':

        # Set the tabix index location
        compressed_dataset.filesystem_idx_location = (
                compressed_dataset.filesystem_location + '.tbi')

        # Make tabix index
        subprocess.check_call([
            TABIX_BINARY, 
            '-p', 'vcf',
            compressed_dataset.get_absolute_location()
        ])
    else:
        # If it's already here, then make sure the index is right
        assert compressed_dataset.filesystem_idx_location == (
                compressed_dataset.filesystem_location + '.tbi'), (
                'Tabix index file location is not correct.')
        assert os.path.exists(
                compressed_dataset.get_absolute_idx_location()), (
                'Tabix index file does not exist on filesystem.')

    return compressed_dataset

def add_bam_file_track(reference_genome, sample_alignment, alignment_type):
    """Update the JBrowse track config file, trackList.json, for this
    ReferenceGenome with a track for the given sample_alignment and alignment_type.
    """
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
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
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
