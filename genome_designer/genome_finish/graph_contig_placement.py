from collections import namedtuple
import os
import re
import subprocess

from django.conf import settings
import networkx as nx
import pysam

from genome_finish.insertion_placement_read_trkg import simple_align_with_bwa_mem
from main.models import Contig
from main.models import Dataset

MAX_DELETION = 500
MAX_DUP = 10
InsertionVertices = namedtuple('InsertionVertices',
        ['exit_ref', 'enter_contig', 'exit_contig', 'enter_ref'])


def graph_contig_placement(contig_list):

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
    simple_align_with_bwa_mem(
            contig_concat,
            get_fasta(ref_genome),
            contig_alignment_bam)

    # Create dictionaries to translate contig uid to its fasta descriptor line
    contig_fasta_descriptor_to_uid = {}
    contig_uid_to_fasta_descriptor = {}
    for contig in contig_list:
        with open(get_fasta(contig), 'r') as fh:
            descriptor = fh.next()
            contig_fasta_descriptor_to_uid[
                    descriptor.strip('>\n')] = contig.uid
            contig_uid_to_fasta_descriptor[
                    contig.uid] = descriptor.strip('>\n')

    # Create graph
    G = nx.DiGraph()

    # Create sequence interval instances for reference and each contig
    ref_intervals = SequenceIntervals(
            ref_genome.uid, ref_genome.num_bases, tag='ref')

    contigs_intervals = {}
    for contig in contig_list:
        m = re.findall('_(\d+)$', contig.label)
        contig_num = m[0]
        tag = 'C' + str(contig_num)
        contigs_intervals[contig.uid] = SequenceIntervals(
                contig.uid, contig.num_bases, tag=tag)

    # Iterate over aligned contig 'reads' in contig alignment to ref bam
    contig_alignmentfile = pysam.AlignmentFile(contig_alignment_bam)
    for read in contig_alignmentfile:

        contig_uid = contig_fasta_descriptor_to_uid[read.qname]
        contig_intervals = contigs_intervals[contig_uid]

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

        # Add chromosome and RC information
        # TODO: Make this consensus based rather than assuming all
        # reads agree
        contig = Contig.objects.get(uid=contig_uid)
        contig.metadata['chromosome'] = contig_alignmentfile.getrname(
                read.reference_id)
        contig.metadata['is_reverse'] = read.is_reverse
        contig.save()

    # Add inter-contig sequence edges
    for contig_intervals in contigs_intervals.values() + [ref_intervals]:
        previous_vertex = contig_intervals.vertices[0]
        for vertex in contig_intervals.vertices[1:]:
            G.add_edge(previous_vertex, vertex)
            previous_vertex = vertex

    placeable_contigs = []
    G.ref_intervals = ref_intervals
    iv_list = novel_seq_ins_walk(G)
    for insertion_vertices in iv_list:
        contig_uid = insertion_vertices.enter_contig.seq_uid
        contig = Contig.objects.get(uid=contig_uid)
        set_contig_placement_params(contig, insertion_vertices)
        placeable_contigs.append(contig)

    return placeable_contigs


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
    for i, exit_ref in enumerate(G.ref_intervals.vertices):
        for enter_contig in contig_neighbors(exit_ref):
            exit_contig = enter_contig
            while exit_contig:
                for enter_ref in ref_neighbors(exit_contig):
                    deletion = enter_ref.pos - exit_ref.pos
                    if -MAX_DUP <= deletion <= MAX_DELETION:
                        iv_list.append(InsertionVertices(
                                exit_ref, enter_contig, exit_contig,
                                enter_ref))
                        break
                exit_contig_list = contig_neighbors(exit_contig)
                exit_contig = (exit_contig_list[0] if len(exit_contig_list)
                        else None)

    return iv_list


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
        cigar = cigar_codes[code]

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
