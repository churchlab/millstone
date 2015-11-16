from collections import defaultdict
import gzip
import os
import re
import subprocess

from django.conf import settings
from settings import JBROWSE_DATA_URL_ROOT
from settings import S3_BUCKET

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from utils.bam_utils import sort_bam_by_coordinate
from utils.bam_utils import index_bam
from utils.import_util import add_dataset_to_entity
from utils.jbrowse_util import prepare_jbrowse_ref_sequence
from utils.jbrowse_util import write_tracklist_json


def maybe_create_reads_to_contig_bam(contig):
    if not contig.dataset_set.filter(
            type=Dataset.TYPE.BWA_ALIGN).exists():
        prepare_jbrowse_ref_sequence(contig)
        align_contig_reads_to_contig(contig)
        add_contig_reads_to_contig_bam_track(contig, Dataset.TYPE.BWA_ALIGN)


def align_contig_reads_to_contig(contig):

    # Get fasta of reads used to make contig
    contig_reads_fasta = os.path.join(
            contig.get_model_data_dir(),
            'extracted_reads.fa')

    # Pull out contig read qnames and put in dictionary contig_reads
    p1 = re.compile('>(\S+)/(\d)')
    contig_reads = defaultdict(list)
    with open(contig_reads_fasta) as fh:
        for line in fh:
            m1 = p1.match(line)
            if m1:
                read_id = m1.group(1)
                read_number = int(m1.group(2))
                contig_reads[read_id].append(read_number)

    # Get source reads fastqs
    sample = contig.experiment_sample_to_alignment.experiment_sample
    source_fq1 = sample.dataset_set.get(
            type=Dataset.TYPE.FASTQ1).get_absolute_location()
    source_fq2_query = sample.dataset_set.filter(
            type=Dataset.TYPE.FASTQ2)
    is_paired_end = source_fq2_query.exists()
    if is_paired_end:
        source_fq2 = source_fq2_query[0].get_absolute_location()

    # Make filenames for contig read fastqs
    output_fq1 = os.path.join(
            contig.get_model_data_dir(),
            'reads.1.fq')
    if is_paired_end:
        output_fq2 = os.path.join(
                contig.get_model_data_dir(),
                'reads.2.fq')

    # Go through source fastqs and write reads in contig_reads to file
    source_fq_list = [source_fq1]
    output_fq_list = [output_fq1]
    if is_paired_end:
        source_fq_list.append(source_fq2)
        output_fq_list.append(output_fq2)

    p1 = re.compile('@(\S+)')
    for input_fq_path, output_fq_path in zip(source_fq_list, output_fq_list):
        if input_fq_path.endswith('.fq'):
            file_like = open(input_fq_path)
        elif input_fq_path.endswith('.gz'):
            file_like = gzip.open(input_fq_path)
        else:
            raise Exception('Compression type not supported')

        with file_like as in_fh, \
             open(output_fq_path, 'w') as out_fh:
            for line in in_fh:
                m1 = p1.match(line)
                if m1:
                    qname = m1.group(1)
                    if qname in contig_reads:
                        out_fh.write(line)
                        out_fh.write(in_fh.next())
                        out_fh.write(in_fh.next())
                        out_fh.write(in_fh.next())

    # Align fastqs to contig fasta
    contig_fasta = contig.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    contig_reads_to_contig_bam = os.path.join(
            contig.get_model_data_dir(),
            'reads_to_contig.bam')
    simple_align_paired_with_bwa_mem(
            output_fq_list,
            contig_fasta,
            contig_reads_to_contig_bam)

    # Coordinate sort and index bam for jbrowse
    coordinate_sorted_bam = (os.path.splitext(contig_reads_to_contig_bam)[0] +
            '.coordinate_sorted.bam')
    sort_bam_by_coordinate(contig_reads_to_contig_bam, coordinate_sorted_bam)
    index_bam(coordinate_sorted_bam)

    # Add the bam file to contig as BWA_ALIGN dataset, overwriting it
    # if it already exists
    dataset_query = contig.dataset_set.filter(
            type=Dataset.TYPE.BWA_ALIGN)
    if dataset_query.count():
        dataset_query[0].delete()

    add_dataset_to_entity(
            contig,
            Dataset.TYPE.BWA_ALIGN,
            Dataset.TYPE.BWA_ALIGN,
            filesystem_location=coordinate_sorted_bam)


def simple_align_paired_with_bwa_mem(reads_fq, reference_fasta,
        output_bam_path):

    # Ensure reference fasta is indexed
    subprocess.call(' '.join([
            '%s/bwa/bwa' % settings.TOOLS_DIR,
            'index',
            reference_fasta
            ]),
            shell=True, executable=settings.BASH_PATH)

    # Align the fastqs to the reference
    align_input_args = [
            '%s/bwa/bwa' % settings.TOOLS_DIR,
            'mem',
            reference_fasta]

    align_input_args.extend(reads_fq)
    align_input_cmd = ' '.join(align_input_args)

    # To skip saving the SAM file to disk directly, pipe output directly to
    # make a BAM file.
    align_input_cmd += (' | ' + settings.SAMTOOLS_BINARY +
            ' view -bS -')

    # Run alignment
    with open(output_bam_path, 'w') as fh:
        subprocess.check_call(
                align_input_cmd, stdout=fh,
                shell=True, executable=settings.BASH_PATH)


def add_contig_reads_to_contig_bam_track(contig, alignment_type):
    """Update the JBrowse track config file, trackList.json, for this
    Contig with a track for the alignment of its reads to itself
    """
    # Get the bam file location from the the Dataset of the contig
    # keyed by the alignment_type.
    bam_dataset = get_dataset_with_type(contig, alignment_type)

    # Figure out the url that JBrowse would use to show the data, e.g.:
    #     /jbrowse/gd_data/projects/58a62c7d/genomes/8dc829ec/align.bam
    # urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
    #         bam_dataset.filesystem_location)

    # NOTE: We should construct bam file urls using project.get_client
    # jbrowse_link() rather than checking S3 flag here.
    if contig.parent_reference_genome.project.is_s3_backed():
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            bam_dataset.filesystem_location.strip("/jbrowse"))
    else:
        urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
            bam_dataset.filesystem_location)

    label = bam_dataset.internal_string(
            contig)

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
    write_tracklist_json(contig, raw_dict_obj, label)

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
    write_tracklist_json(contig, coverage_raw_dict_obj,
            snp_coverage_label)


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
