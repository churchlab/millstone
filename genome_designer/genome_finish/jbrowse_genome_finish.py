from collections import defaultdict
import os
import re

from settings import JBROWSE_DATA_URL_ROOT
from settings import S3_BUCKET

from main.model_utils import get_dataset_with_type
from main.models import Dataset

from utils.jbrowse_util import write_tracklist_json


def align_contig_reads_to_contig(contig):

    # Get reads
    contig_reads_fasta = os.path.join(
            contig.get_model_data_dir(),
            'extracted_reads.fa')

    p1 = re.compile('>(\S+)/(\d)')
    contig_reads = defaultdict(list)
    with open(contig_reads_fasta) as fh:
        for line in fh:
            m1 = p1.match(line)
            if m1:
                read_id = m1.group(1)
                read_number = int(m1.group(2))
                contig_reads[read_id].append(read_number)

    sample = contig.experiment_sample_to_alignment.experiment_sample
    source_fq1 = sample.dataset_set.get(
            Dataset.TYPE.FASTQ1).get_absolute_location()
    source_fq2 = sample.dataset_set.get(
            Dataset.TYPE.FASTQ2).get_absolute_location()

    output_fq1 = os.path.join(
            contig.get_model_data_dir(),
            'reads.1.fq')
    output_fq2 = os.path.join(
            contig.get_model_data_dir(),
            'reads.2.fq')

    # Go through fastqs and write reads in ROI to file
    p1 = re.compile('@(\S+)')
    for input_fq_path, output_fq_path in [
            (source_fq1, output_fq1), (source_fq2, output_fq2)]:
        counter = 0
        if desired_coverage:
            iterator = iter(include_read)
        with open(input_fq_path, 'r') as in_fh, open(output_fq_path, 'w') as out_fh:
            for line in in_fh:
                m1 = p1.match(line)
                if m1:
                    qname = m1.group(1)
                    if qname in qnames_in_region:
                        if desired_coverage:
                            if iterator.next():
                                out_fh.write(line)
                                out_fh.write(in_fh.next())
                                out_fh.write(in_fh.next())
                                out_fh.write(in_fh.next())
                        else:
                            out_fh.write(line)
                            out_fh.write(in_fh.next())
                            out_fh.write(in_fh.next())
                            out_fh.write(in_fh.next())

def write_bam_to_fastqs(reads, fastq_path):
    """Writes the aligned portion of each read into a fastq
    """

    query_alignment_seqrecords = []
    for read in reads:
        query_alignment_seqrecords.append(SeqRecord(
                Seq(read.query_alignment_sequence, IUPAC.ambiguous_dna),
                letter_annotations={
                        'phred_quality': read.query_alignment_qualities},
                id=read.query_name,
                description=''))

    with open(fastq_path, 'w') as fastq_handle:
        SeqIO.write(query_alignment_seqrecords, fastq_handle, 'fastq')


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
                'chunkSizeLimit': 10000000, # double the default chunk size
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
