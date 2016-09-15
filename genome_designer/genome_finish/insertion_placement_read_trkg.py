from collections import defaultdict
import os
import pickle
import re
import subprocess

from django.conf import settings
import pysam

from main.models import Contig
from main.models import Dataset
from pipeline.read_alignment_util import has_bwa_index
from utils.bam_utils import index_bam
from utils.bam_utils import sort_bam_by_coordinate
from utils.import_util import add_dataset_to_entity


def mapped_mates_of_unmapped_reads(contig):
    unmapped_contig_reads = extract_contig_reads(
            contig, read_category='unmapped')
    print len(unmapped_contig_reads), 'unmapped reads in contig'

    original_align = contig.experiment_sample_to_alignment.dataset_set.get(
            type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    original_alignmentfile = pysam.AlignmentFile(original_align)
    found_mates = []
    for read in unmapped_contig_reads:
        if not read.mate_is_unmapped:
            mate = original_alignmentfile.mate(read)
            found_mates.append(mate)
    original_alignmentfile.close()

    print len(found_mates), 'mapped mates found'
    return found_mates


def extract_contig_reads(contig, read_category='all'):
    '''
    Use velvet tools to extract reads from the contig, and then
    pull them out of the original reference genome alignment.
    It returns a list of the reads as pysam objects and also
    puts them as a 'extracted_reads.fa' in the contig data dir.
    '''

    READ_CATEGORY_TO_FILENAME_DICT = {
        'without_mates': 'bwa_align.SV_indicants_no_dups.bam',
        'clipped': 'bwa_align.clipped.bam',
        'split': 'bwa_align.split.bam',
        'unmapped': 'bwa_align.unmapped.bam'
    }

    def _read_category_to_filename(read_category):
        if read_category in READ_CATEGORY_TO_FILENAME_DICT:
            return READ_CATEGORY_TO_FILENAME_DICT[read_category]
        elif read_category == 'all':

            assembly_metadata_file = os.path.join(
                    contig.metadata['assembly_dir'],
                    'metadata.txt')
            with open(assembly_metadata_file) as fh:
                assembly_metadata_obj = pickle.load(fh)
            return assembly_metadata_obj['sv_indicants_bam']
        elif read_category == 'mates_of_unmapped':
            return mapped_mates_of_unmapped_reads(contig)
        else:
            raise Exception('read category not recognized')

    extract_contig_reads_executable = os.path.join(
            settings.TOOLS_DIR,
            'velvet/extractContigReads.pl')

    assembly_dir = contig.metadata['assembly_dir']

    contig_node_number = contig.metadata['node_number']
    cmd = [extract_contig_reads_executable, str(contig_node_number),
           assembly_dir]
    cmd = ' '.join(cmd)

    contig_reads_fasta = os.path.join(
            contig.get_model_data_dir(),
            'extracted_reads.fa')
    if not os.path.exists(contig_reads_fasta):
        with open(contig_reads_fasta, 'w') as fh:
            subprocess.call(cmd, shell=True, stdout=fh)

    p1 = re.compile('>(\S+)/(\d)')
    contig_reads = defaultdict(list)
    with open(contig_reads_fasta) as fh:
        for line in fh:
            m1 = p1.match(line)
            if m1:
                read_id = m1.group(1)
                read_number = int(m1.group(2))
                contig_reads[read_id].append(read_number)

    sv_indicant_reads_path = os.path.join(
            contig.experiment_sample_to_alignment.get_model_data_dir(),
            _read_category_to_filename(read_category))

    sam_file = pysam.AlignmentFile(sv_indicant_reads_path)
    sv_indicant_reads_in_contig = []
    for read in sam_file:
        if read.is_read1:
            read_number = 1
        elif read.is_read2:
            read_number = 2
        else:
            raise Exception('Read is neither read1 nor read2')

        contig_read_numbers = contig_reads.get(read.query_name, [])
        if read_number in contig_read_numbers:
            sv_indicant_reads_in_contig.append(read)

    # HACK: Set chromosome here while sam file is open
    # so AlignmentFile.getrname(tid) can be called
    ref_id_to_count = {}
    mapped_count = 0
    for read in sv_indicant_reads_in_contig:
        if not read.is_unmapped:
            mapped_count += 1
            if read.reference_id not in ref_id_to_count:
                ref_id_to_count[read.reference_id] = 1
            else:
                ref_id_to_count[read.reference_id] += 1

    if mapped_count:
        tid_count_sorted = sorted(
                ref_id_to_count.items(), key=lambda x: x[1], reverse=True)

        mode_chrom_tid = tid_count_sorted[0][0]
        mode_chrom_percentage = (tid_count_sorted[0][1] /
                float(mapped_count))

        # Set field
        if mode_chrom_percentage > 0.8:
            contig_seqrecord_id = sam_file.getrname(mode_chrom_tid)
            contig.metadata['chromosome'] = contig_seqrecord_id
            contig.save()

    sam_file.close()
    return sv_indicant_reads_in_contig


def make_contig_reads_dataset(contig, sv_indicant_reads_in_contig):
    '''
    Using the contig reads generated by extract_contig_reads(),
    generate a bam file, index and sort it.
    '''
    # Get bam filename
    extracted_reads_bam_file = os.path.join(
            contig.get_model_data_dir(),
            'sv_indicants.bam')

    bwa_align_bam = contig.experiment_sample_to_alignment.dataset_set.get(
            type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    sam_file = pysam.AlignmentFile(bwa_align_bam)

    # Write extracted reads into bam file
    extracted_reads_alignment_file = pysam.AlignmentFile(
            extracted_reads_bam_file, "wb", template=sam_file)
    sam_file.close()

    for read in sv_indicant_reads_in_contig:
        extracted_reads_alignment_file.write(read)

    extracted_reads_alignment_file.close()

    coordinate_sorted_bam = (os.path.splitext(extracted_reads_bam_file)[0] +
            '.coordinate_sorted.bam')
    sort_bam_by_coordinate(extracted_reads_bam_file, coordinate_sorted_bam)
    index_bam(coordinate_sorted_bam)

    # Add the bam file to contig as BWA_SV_INDICANTS dataset, overwriting it
    # if it already exists
    dataset_query = contig.dataset_set.filter(
            type=Dataset.TYPE.BWA_SV_INDICANTS)
    if dataset_query.count():
        dataset_query[0].delete()

    add_dataset_to_entity(contig,
            Dataset.TYPE.BWA_SV_INDICANTS,
            Dataset.TYPE.BWA_SV_INDICANTS,
            filesystem_location=coordinate_sorted_bam)


def extract_left_and_right_clipped_read_dicts(sv_indicant_reads_in_contig,
        clipping_threshold=0):

    SOFT_CLIP = 4
    HARD_CLIP = 5
    CLIP = [SOFT_CLIP, HARD_CLIP]

    # Separate left and right clipped reads
    left_clipped = defaultdict(list)
    right_clipped = defaultdict(list)
    for read in sv_indicant_reads_in_contig:
        if read.cigartuples is not None:
            left_clipping = (read.cigartuples[0][1]
                    if read.cigartuples[0][0] in CLIP else 0)
            right_clipping = (read.cigartuples[-1][1]
                    if read.cigartuples[-1][0] in CLIP else 0)
            if max(left_clipping, right_clipping) > clipping_threshold:
                is_left_clipped = left_clipping > right_clipping
                is_right_clipped = right_clipping > left_clipping
                if is_left_clipped:
                    left_clipped[read.reference_start].append(read)
                elif is_right_clipped:
                    right_clipped[read.reference_end].append(read)

    return {
        'left_clipped': left_clipped,
        'right_clipped': right_clipped
    }


def simple_align_with_bwa_mem(reads_fq,
        reference_fasta, output_bam_path, bwa_arg_list=[]):

    # Assert reference fasta is indexed
    assert has_bwa_index(reference_fasta)

    # Align clipped query alignment fastq to contig
    align_input_args = ['%s/bwa/bwa' % settings.TOOLS_DIR, 'mem']

    # add mem args
    align_input_args.extend(bwa_arg_list)
    align_input_args.extend([reference_fasta, reads_fq])

    align_cmd = ' '.join(align_input_args)

    # Bwa mem calls reads clipped slightly at the end of the genome
    # as unmapped, so filter these out with -F 0x004
    # To skip saving the SAM file to disk directly, pipe output directly to
    # make a BAM file.
    align_cmd += (' | ' + settings.SAMTOOLS_BINARY +
            ' view -F 0x004 -bS -')

    # Run alignment
    with open(output_bam_path, 'w') as fh:
        subprocess.check_call(
                align_cmd, stdout=fh,
                shell=True, executable=settings.BASH_PATH)


def make_contig_reads_to_ref_alignments(
        contig_uid, add_jbrowse_track=False, overwrite=False):
    """NOTE: add_jbrowse_track not currently supported.
    """
    if add_jbrowse_track:
        raise NotImplementedError("JBrowse tracks disabled.")

    contig = Contig.objects.get(uid=contig_uid)

    dataset_query = contig.dataset_set.filter(
            type=Dataset.TYPE.BWA_SV_INDICANTS)

    if overwrite or not dataset_query.count():
        # Get the reads aligned to ref that assembled the contig
        contig_reads = extract_contig_reads(contig, 'all')

        if len(contig_reads) == 0:
            raise Exception(
                    'No reads were extracted from contig ' + contig.label)

        # Add contig reads to contig dataset set as dataset
        # type BWA_SV_INDICANTS
        make_contig_reads_dataset(contig, contig_reads)

    # TODO(dbgoodman): Maybe fix.
    # if add_jbrowse_track:
    #     add_contig_reads_bam_track(contig, Dataset.TYPE.BWA_SV_INDICANTS)
