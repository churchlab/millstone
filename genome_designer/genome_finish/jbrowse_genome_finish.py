import os

from settings import JBROWSE_DATA_URL_ROOT
from settings import S3_BUCKET

from main.model_utils import get_dataset_with_type
from utils.jbrowse_util import write_tracklist_json


def add_contig_reads_bam_track(contig, alignment_type):
    """Update the JBrowse track config file, trackList.json, for this
    ReferenceGenome with a track for the given sample_alignment and
    alignment_type.
    """
    # Get the bam file location from the the Dataset of the genome
    # keyed by the alignment_type.
    bam_dataset = get_dataset_with_type(contig, alignment_type)

    # assert bam_dataset.is_indexed()

    # Figure out the url that JBrowse would use to show the data, e.g.:
    #     /jbrowse/gd_data/projects/58a62c7d/genomes/8dc829ec/align.bam
    # urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
    #         bam_dataset.filesystem_location)

    # NOTE: We should construct bam file urls using project.get_client
    # jbrowse_link() rather than checking S3 flag here.
    reference_genome = contig.parent_reference_genome
    if reference_genome.project.is_s3_backed():
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            bam_dataset.filesystem_location.strip("/jbrowse"))
    else:
        urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
            bam_dataset.filesystem_location)

    # doing label as ES_AG because SA isn't currently used in the variant view
    label = bam_dataset.internal_string(contig)

    key = bam_dataset.external_string(contig)

    # Build the JSON object.
    raw_dict_obj = {
        'tracks': [
            {
                'storeClass': 'JBrowse/Store/SeqFeature/BAM',
                'urlTemplate': urlTemplate,
                'label': label,
                'type': 'JBrowse/View/Track/Alignments2',
                'chuckSizeLimit': 10000000, # double the default chunk size
                'key': key,
                'category': 'Contig BAM Tracks',
                'style': {
                    'className': 'alignment',
                    'arrowheadClass': 'arrowhead',
                    'labelScale': 100
                }
            }
        ]}
    write_tracklist_json(reference_genome, raw_dict_obj, label)

    # Also add a snp coverage track.
    snp_coverage_label = bam_dataset.internal_string(
            contig) + '_COVERAGE'

    snp_coverage_key = key + ' Coverage'
    coverage_raw_dict_obj = {
        'tracks': [
            {
                'storeClass': 'JBrowse/Store/SeqFeature/BAM',
                'urlTemplate': urlTemplate,
                'label': snp_coverage_label,
                'type': 'JBrowse/View/Track/SNPCoverage',
                'category': 'Contig Coverage Tracks',
                'key': snp_coverage_key
            }
        ]}
    write_tracklist_json(reference_genome, coverage_raw_dict_obj,
            snp_coverage_label)
