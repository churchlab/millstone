from collections import namedtuple
import os
import subprocess

from Bio import SeqIO
from django.conf import settings
import networkx as nx
import pysam

from genome_finish.contig_display_utils import Junction
from genome_finish.millstone_de_novo_fns import get_coverage_stats
from genome_finish.insertion_placement_read_trkg import extract_contig_reads
from genome_finish.insertion_placement_read_trkg import make_contig_reads_dataset
from genome_finish.insertion_placement_read_trkg import simple_align_with_bwa_mem
from genome_finish.jbrowse_genome_finish import add_contig_reads_bam_track
from main.models import Contig
from main.models import Dataset
from main.models import ReferenceGenome

MAX_DELETION = 50000
MAX_DUP = 10
MAX_REF_SELF_HOMOLOGY = 200
InsertionVertices = namedtuple('InsertionVertices',
        ['exit_ref', 'enter_contig', 'exit_contig', 'enter_ref'])


def graph_contig_placement(contig_list, skip_extracted_read_alignment,
        use_alignment_reads=True):

    if not skip_extracted_read_alignment:
        # Make a bam track on the reference for each contig that shows only the
        # reads that assembled the contig and their mates
        for contig in contig_list:
            make_contig_reads_to_ref_alignments(contig)

    # Make Assembly dir
    assembly_dir = contig_list[0].metadata['assembly_dir']
    contig_alignment_dir = os.path.join(assembly_dir, 'contig_alignment')
    os.mkdir(contig_alignment_dir)

    # Concatenate contig fastas for alignment
    contig_fastas = [get_fasta(c) for c in contig_list]
    contig_concat = os.path.join(contig_alignment_dir, 'contig_concat.fa')
    with open(contig_concat, 'w') as fh:
        subprocess.check_call(' '.join(['cat'] + contig_fastas),
                shell=True, executable=settings.BASH_PATH, stdout=fh)

    # Align concatenated contig fastas to reference
    ref_genome = contig_list[0].parent_reference_genome
    contig_alignment_bam = os.path.join(
            contig_alignment_dir, 'contig_alignment.bam')
    print 'Aligning contigs to reference'
    simple_align_with_bwa_mem(
            contig_concat,
            get_fasta(ref_genome),
            contig_alignment_bam)

    # Create graph
    G = nx.DiGraph()

    # Create sequence interval instances for reference and each contig
    ref_intervals = SequenceIntervals(
            ref_genome.uid, ref_genome.num_bases, tag='ref')

    G.ref_intervals = ref_intervals

    add_alignment_to_graph(G, contig_alignment_bam)
    detect_strand_chromosome_junctions(contig_list, contig_alignment_bam)

    # Create dictionaries to translate contig uid to its fasta descriptor line
    contig_qname_to_uid = {}
    for contig in contig_list:
        with open(get_fasta(contig), 'r') as fh:
            descriptor = fh.next()
            contig_qname_to_uid[
                    descriptor.strip('>\n')] = contig.uid

    placeable_contigs = []
    iv_list = novel_seq_ins_walk(G)
    sample_alignment = contig_list[0].experiment_sample_to_alignment
    if use_alignment_reads:
        coverage_stats = get_coverage_stats(sample_alignment)
        sample_alignment_bam = sample_alignment.dataset_set.get(
            type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    for insertion_vertices in iv_list:

        contig_qname = insertion_vertices.enter_contig.seq_uid
        contig_uid = contig_qname_to_uid[contig_qname]
        contig = Contig.objects.get(uid=contig_uid)
        set_contig_placement_params(contig, insertion_vertices)

        if use_alignment_reads:
            # Filter out deletions of good coverage regions
            deletion_length = (insertion_vertices.enter_ref.pos -
                    insertion_vertices.exit_ref.pos)

            if deletion_length > 0:
                deletion_cov = avg_coverage(
                        sample_alignment_bam,
                        contig.metadata['chromosome'],
                        insertion_vertices.exit_ref.pos,
                        insertion_vertices.enter_ref.pos)

            chrom_cov_stats = coverage_stats[contig.metadata['chromosome']]
            chrom_cov_mean = chrom_cov_stats['mean']
            chrom_cov_std = chrom_cov_stats['std']
            if deletion_length <= 0 or (
                    deletion_cov < chrom_cov_mean - chrom_cov_std):
                placeable_contigs.append(contig)
        else:
            placeable_contigs.append(contig)

    # Perform translocation walk
    if ref_genome.num_chromosomes == 1:

        trans_iv_pairs = translocation_walk(G)
        var_dict_list = [parse_path_into_ref_alt(iv_pair, contig_qname_to_uid)
                for iv_pair in trans_iv_pairs]
    else:
        print 'Translocation walk not implemented for multi-chromosomal refs'
        var_dict_list = []

    return placeable_contigs, var_dict_list


def avg_coverage(bam, chromosome, start, end):

    alignment_file = pysam.AlignmentFile(bam)
    coverage_atcg = alignment_file.count_coverage(
            str(chromosome), start, end)

    avg_coverage = sum(max(atcg) for atcg in zip(*[coverage_atcg[i] for i in
            range(4)])) / (end - start)

    return avg_coverage


def add_alignment_to_graph(G, contig_alignment_bam):

    ref_intervals = G.ref_intervals

    # Iterate over aligned contig 'reads' in contig alignment to ref bam
    contigs_intervals = {}
    contig_alignmentfile = pysam.AlignmentFile(contig_alignment_bam)
    for read in contig_alignmentfile:
        if read.is_secondary or read.is_supplementary:
            continue

        assert read.qname not in contigs_intervals
        assert 'H' not in read.cigarstring

        contigs_intervals[read.qname] = SequenceIntervals(
                read.qname, read.query_length)

    # Iterate over aligned contig 'reads' in contig alignment to ref bam
    contig_alignmentfile = pysam.AlignmentFile(contig_alignment_bam)
    for read in contig_alignmentfile:

        contig_intervals = contigs_intervals[read.qname]

        match_regions = get_match_regions(read)
        for match_region in match_regions:
            print match_region
            print read.query_length
            if match_region.read_start != 0:
                # Insert break into sequence intervals for ref and contig
                ref_vert = ref_intervals.insert_vertex(match_region.ref_start)
                contig_vert = contig_intervals.insert_vertex(
                        match_region.read_start)

                # Add an edge to G for the junction
                G.add_edge(contig_vert, ref_vert, weight=match_region.length)

            if match_region.read_end != contig_intervals.length:
                # Insert break into sequence intervals for ref and contig
                ref_vert = ref_intervals.insert_vertex(match_region.ref_end)
                contig_vert = contig_intervals.insert_vertex(
                        match_region.read_end)

                # Add an edge to G for the junction
                G.add_edge(ref_vert, contig_vert, weight=match_region.length)

    # Add inter-contig sequence edges
    for contig_intervals in contigs_intervals.values() + [ref_intervals]:
        previous_vertex = contig_intervals.vertices[0]
        for vertex in contig_intervals.vertices[1:]:
            G.add_edge(previous_vertex, vertex)
            previous_vertex = vertex

    # Add back edges in contigs
    for contig_intervals in contigs_intervals.values():
        previous_vertex = contig_intervals.vertices[0]
        for vertex in contig_intervals.vertices[1:]:
            G.add_edge(vertex, previous_vertex)
            previous_vertex = vertex

    # G.ref_intervals = ref_intervals
    G.contig_intervals_list = contigs_intervals


def detect_strand_chromosome_junctions(contig_list, contig_alignment_bam):
    # Create dictionaries to translate contig uid to its fasta descriptor line
    contig_qname_to_uid = {}
    for contig in contig_list:
        with open(get_fasta(contig), 'r') as fh:
            descriptor = fh.next()
            contig_qname_to_uid[
                    descriptor.strip('>\n')] = contig.uid

    # Iterate over aligned contig 'reads' in contig alignment to ref bam
    contig_alignmentfile = pysam.AlignmentFile(contig_alignment_bam)
    for read in contig_alignmentfile:

        contig_uid = contig_qname_to_uid[read.qname]

        match_regions = get_match_regions(read)

        # Add chromosome and RC information
        # TODO: Make this consensus based rather than assuming all
        # reads agree
        contig = Contig.objects.get(uid=contig_uid)
        contig.metadata['chromosome'] = contig_alignmentfile.getrname(
                read.reference_id)
        contig.metadata['is_reverse'] = read.is_reverse
        set_contig_junctions(contig, match_regions)
        contig.save()


def set_contig_junctions(contig, match_region_list):
    left_junctions = []
    right_junctions = []
    for i, match_region in enumerate(match_region_list):
        if match_region.read_start != 0:

            prev_match = match_region_list[i-1] if i != 0 else None
            if prev_match:
                new_seq_length = match_region.read_start - prev_match.read_end
            else:
                new_seq_length = match_region.read_start

            right_junctions.append(Junction(
                    match_region.ref_start, match_region.length,
                    match_region.read_start, new_seq_length))

        if match_region.read_end != contig.num_bases:

            next_match = (match_region_list[i+1] if
                    i != len(match_region_list) - 1 else None)
            if next_match:
                new_seq_length = match_region.read_end - next_match.read_start
            else:
                new_seq_length = contig.num_bases - match_region.read_end

            left_junctions.append(Junction(
                    match_region.ref_end, match_region.length,
                    match_region.read_end, new_seq_length))

    for key, data in [('left_junctions', left_junctions),
            ('right_junctions', right_junctions)]:

        if isinstance(contig.metadata.get(key), list):
            contig.metadata[key].extend(data)
        else:
            contig.metadata[key] = data
    contig.save()


def set_contig_placement_params(contig, insertion_vertices):
    contig.metadata['contig_insertion_endpoints'] = (
            insertion_vertices.enter_contig.pos,
            insertion_vertices.exit_contig.pos)
    contig.metadata['reference_insertion_endpoints'] = (
            insertion_vertices.exit_ref.pos,
            insertion_vertices.enter_ref.pos)
    contig.save()


def novel_seq_ins_walk(G):
    ref_seq_uid = G.ref_intervals.vertices[0].seq_uid

    def ref_neighbors(vert):
        return [v for v in G.neighbors(vert) if v.seq_uid == ref_seq_uid]

    def contig_neighbors(vert):
        return [v for v in G.neighbors(vert) if v.seq_uid != ref_seq_uid]

    iv_list = []
    for exit_ref in G.ref_intervals.vertices:
        for enter_contig in contig_neighbors(exit_ref):
            queue = [enter_contig]
            visited = []
            while queue:
                exit_contig = queue.pop()
                for enter_ref in ref_neighbors(exit_contig):
                    deletion = enter_ref.pos - exit_ref.pos
                    ref_self_homology = enter_contig.pos - exit_contig.pos
                    if (-MAX_DUP < deletion < MAX_DELETION and
                            ref_self_homology < MAX_REF_SELF_HOMOLOGY):
                        iv_list.append(InsertionVertices(
                                exit_ref, enter_contig, exit_contig,
                                enter_ref))
                        break
                visited.append(exit_contig)
                queue.extend([n for n in contig_neighbors(exit_contig)
                        if n not in visited])

    return iv_list


def translocation_walk(G):
    ref_seq_uid = G.ref_intervals.vertices[0].seq_uid

    def ref_neighbors(vert):
        return [v for v in G.neighbors(vert) if v.seq_uid == ref_seq_uid]

    def contig_neighbors(vert):
        return [v for v in G.neighbors(vert) if v.seq_uid != ref_seq_uid]

    forward_edges = []
    back_edges = []

    for exit_ref in G.ref_intervals.vertices:
        for enter_contig in contig_neighbors(exit_ref):
            queue = [enter_contig]
            visited = []
            while queue:
                exit_contig = queue.pop()
                for enter_ref in ref_neighbors(exit_contig):

                    iv = InsertionVertices(
                            exit_ref, enter_contig, exit_contig,
                            enter_ref)

                    if exit_ref.pos < enter_ref.pos:
                        forward_edges.append(iv)
                    else:
                        back_edges.append(iv)

                visited.append(exit_contig)
                queue.extend([n for n in contig_neighbors(exit_contig)
                        if n not in visited])

    sorted_by_exit_ref = sorted(forward_edges + back_edges,
            key=lambda x: x.exit_ref.pos)

    sorted_by_enter_ref = sorted(forward_edges + back_edges,
            key=lambda x: x.enter_ref.pos)

    MAX_TRANS_LENGTH = 20000
    MIN_TRANS_LENGTH = 20

    iv_pairs = []
    i = 0
    for enter_iv in sorted_by_exit_ref:
        while (sorted_by_enter_ref[i].enter_ref.pos < enter_iv.exit_ref.pos and
                i < len(sorted_by_enter_ref) -1 ):
            i += 1

        j = i
        exit_iv = sorted_by_enter_ref[j]
        deletion = exit_iv.enter_ref.pos - enter_iv.exit_ref.pos
        while -MAX_DUP < deletion < MAX_DELETION:

            trans_length = exit_iv.exit_ref.pos - enter_iv.enter_ref.pos
            if MIN_TRANS_LENGTH < trans_length < MAX_TRANS_LENGTH:
                iv_pairs.append((enter_iv, exit_iv))

            if j == len(sorted_by_enter_ref) - 1:
                break

            j += 1
            exit_iv = sorted_by_enter_ref[j]
            deletion = exit_iv.enter_ref.pos - enter_iv.exit_ref.pos

    return iv_pairs


class SequenceVertex:

    def __init__(self, seq_uid, pos, parent):
        self.parent = parent
        self.seq_uid = seq_uid
        self.pos = pos
        self.uid = '_'.join([self.seq_uid, str(self.pos)])

    def __repr__(self):
        return self.parent.tag + '_' + str(self.pos)

    def __hash__(self):
        return hash(self.uid)

    def __cmp__(self, rhs):
        if not isinstance(rhs, SequenceVertex):
            return TypeError
        if rhs.seq_uid != self.seq_uid:
            return ValueError

        if self.pos < rhs.pos:
            return -1
        if self.pos > rhs.pos:
            return 1
        return 0


class SequenceIntervals:

    def __init__(self, seq_uid, length, tag=None):
        self.seq_uid = seq_uid
        self.tag = tag if tag is not None else seq_uid
        self.length = length
        self.vertices = [SequenceVertex(seq_uid, 0, self),
                         SequenceVertex(seq_uid, length, self)]

    def insert_vertex(self, pos):
        for i, vertex in enumerate(self.vertices):
            if vertex.pos < pos:
                continue
            elif vertex.pos == pos:
                return vertex
            elif vertex.pos > pos:
                new_vertex = SequenceVertex(self.seq_uid, pos, self)
                self.vertices.insert(i, new_vertex)
                return new_vertex


def get_fasta(has_fasta):
    return has_fasta.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()


def get_match_regions(read):

    cigar_codes = {
        0: 'M',
        1: 'I',
        2: 'D',
        4: 'C',
        5: 'C'
    }

    ref_pos = read.reference_start
    read_pos = 0
    matching_length = 0
    matching = False
    MatchRegion = namedtuple('MatchRegion',
            ['ref_start', 'ref_end', 'read_start', 'read_end', 'length'])
    match_regions = []
    print read.cigarstring
    for code, length in read.cigartuples:
        cigar = cigar_codes.get(code, None)
        if cigar is None:
            print 'Cigar code:', code, 'not recognized in read', read
            continue

        if cigar == 'M' and not matching:
            match_ref_start = ref_pos
            match_read_start = read_pos
            matching = True

        elif cigar == 'C' and matching:
            match_regions.append(MatchRegion(
                        match_ref_start, ref_pos,
                        match_read_start, read_pos,
                        matching_length))
            matching = False
            matching_length = 0

        if cigar == 'C':
            read_pos += length
        elif cigar == 'M':
            matching_length += length
            read_pos += length
            ref_pos += length
        elif cigar == 'I':
            read_pos += length
        elif cigar == 'D':
            ref_pos += length

    if matching:
        match_regions.append(MatchRegion(
                    match_ref_start, ref_pos,
                    match_read_start, read_pos,
                    matching_length))

    return match_regions


def make_contig_reads_to_ref_alignments(contig):

    # Get the reads aligned to ref that assembled the contig
    contig_reads = extract_contig_reads(contig, 'all')

    if len(contig_reads) == 0:
        raise Exception('No reads were extracted from contig ' + contig.label)

    # Add contig reads to contig dataset set as dataset type BWA_SV_INDICANTS
    make_contig_reads_dataset(contig, contig_reads)

    # Add bam track
    add_contig_reads_bam_track(contig, Dataset.TYPE.BWA_SV_INDICANTS)

def parse_path_into_ref_alt(path_list, contig_qname_to_uid):

    ref_uid = path_list[0].exit_ref.seq_uid
    ref_genome = ReferenceGenome.objects.get(uid=ref_uid)
    ref_fasta = get_fasta(ref_genome)
    with open(ref_fasta) as fh:
        ref_seqrecord = SeqIO.parse(fh, 'fasta').next()
        ref_seq = str(ref_seqrecord.seq)
        ref_chromosome = ref_seqrecord.id

    def _seq_str(enter_vert, exit_vert):
        if enter_vert.seq_uid == ref_uid:
            return ref_seq[enter_vert.pos: exit_vert.pos]

        contig_qname = enter_vert.seq_uid
        contig_uid = contig_qname_to_uid[contig_qname]
        contig = Contig.objects.get(uid=contig_uid)
        contig_fasta = get_fasta(contig)
        with open(contig_fasta) as fh:
            contig_seqrecord = SeqIO.parse(fh, 'fasta').next()

        # Determine whether contig is reverse complement relative to reference
        is_reverse = contig.metadata.get('is_reverse', False)

        # Extract cassette sequence from contig
        if is_reverse:
            return str(contig_seqrecord.seq.reverse_complement()[
                    enter_vert.pos:exit_vert.pos])
        else:
            return str(contig_seqrecord.seq[
                    enter_vert.pos:exit_vert.pos])

    path_list_concat = reduce(lambda x, y: x + y, path_list)
    def _seq_interval_iter():
        i = 1
        while i <= len(path_list_concat) - 3:
            yield path_list_concat[i:i+2]
            i += 2

    seq_list = []
    for enter_vert, exit_vert in _seq_interval_iter():

        seq_len = exit_vert.pos - enter_vert.pos

        if seq_len < 0:
            seq_list.append(seq_len)
        elif seq_len > 0:
            seq_list.append(_seq_str(enter_vert, exit_vert))
        else:
            seq_list.append('')

    ref_start = path_list[0].exit_ref.pos
    ref_end = path_list[-1].enter_ref.pos
    is_first = True
    alt_seq = ''
    for seq in seq_list:
        if is_first:
            if isinstance(seq, int):
                ref_start += seq
            else:
                alt_seq += seq
                is_first = False
        else:
            if isinstance(seq, int):
                alt_seq = alt_seq[:seq]
            else:
                alt_seq += seq

    ref_seq = ref_seq[ref_start: ref_end]

    return {
        'chromosome': ref_chromosome,
        'pos': ref_start,
        'ref_seq': ref_seq,
        'alt_seq': alt_seq
    }
