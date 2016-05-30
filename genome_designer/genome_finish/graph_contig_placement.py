from collections import namedtuple
import os
import subprocess
import sys

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
from main.models import ExperimentSampleToAlignment
from pipeline.read_alignment_util import ensure_bwa_index
from utils.import_util import add_dataset_to_entity

MAX_DELETION = 100000
MAX_DUP = 10
MAX_REF_SELF_HOMOLOGY = 200
MAX_TRANS_LENGTH = 20000
MIN_TRANS_LENGTH = 20
MAX_TRANS_DELETION = 20000
InsertionVertices = namedtuple('InsertionVertices',
        ['exit_ref', 'enter_contig', 'exit_contig', 'enter_ref'])


def get_genbank_features_with_type(genbank_path, feature_type):

    chrom_intervals = {}
    with open(genbank_path, 'r') as fh:
        for seq_record in SeqIO.parse(fh, 'genbank'):
            interval_list = []
            for f in seq_record.features:
                if f.type == feature_type:
                    interval_list.append((f, f.extract(seq_record.seq)))

            chrom_intervals[seq_record.id] = interval_list

    return chrom_intervals


def write_me_features_multifasta(genbank_path, output_fasta_path,
        also_inversions=False):
    """Write mobile elements to file
    """
    me_features = get_genbank_features_with_type(
            genbank_path, 'mobile_element')
    me_sequences = {}
    for chrom, feature_seq_tuples in me_features.items():
        for feature, seq in feature_seq_tuples:
            if 'mobile_element_type' not in feature.qualifiers:
                continue
            me_type = feature.qualifiers['mobile_element_type'][0]
            me_type = me_type.replace(' ', '_')
            if me_type not in me_sequences:
                me_sequences[me_type] = seq

    with open(output_fasta_path, 'w') as fh:
        for me_type, seq in me_sequences.items():
            fh.write('>ME_%s\n' % me_type)
            fh.write(str(seq))
            fh.write('\n')

        if also_inversions:
            for me_type, seq in me_sequences.items():
                fh.write('>ME_REVERSE_%s\n' % me_type)
                fh.write(str(seq)[::-1])
                fh.write('\n')


def graph_contig_placement(contig_list, skip_extracted_read_alignment,
        use_alignment_reads=True):
    """Align contigs passed in contig_list to the reference and to any
    annotated mobile elements in the reference genbank and use the alignment
    to build a sequence graph.  The sequence graph is then used by graph
    walking algorithms that call structural variants from paths in the graph.

    Args:
        contig_list: list of Contig objects
        skip_extracted_read_alignment: if False, extract the reads that
                assembled each contig and make them a bam track on the
                reference
        use_alignment_reads: if True, filter contig placements that would
                delete regions of moderate coverage

    Returns:
         placeable_contigs: Contig objects with metadata fields holding
                their reference placement parameters
         var_dict_list: list of dictionary representations of translocation
                variants with keys: chromosome, pos, ref_seq, alt_seq
         me_var_dict_list: list of dictionary representations of mobile
                element translocation variants with keys: chromosome, pos,
                ref_seq, alt_seq, MEINFO
    """

    sample_alignment = contig_list[0].experiment_sample_to_alignment
    sample_alignment.data['assembly_status'] = (
            ExperimentSampleToAlignment.ASSEMBLY_STATUS.BUILDING_SEQUENCE_GRAPH
    )
    sample_alignment.save()

    if not skip_extracted_read_alignment:
        # Make a bam track on the reference for each contig that shows only the
        # reads that assembled the contig and their mates
        for contig in contig_list:
            make_contig_reads_to_ref_alignments(contig)

    # Make Assembly dir
    assembly_dir = contig_list[0].metadata['assembly_dir']

    contig_alignment_base_dir = os.path.join(
            assembly_dir, 'contig_alignment')
    if not os.path.exists(contig_alignment_base_dir):
        os.mkdir(contig_alignment_base_dir)

    # Make a subdirectory within contig_alignment_base_dir with the
    # next available integer that is not already a subdir as its name
    contig_alignment_dir = (
            os.path.join(contig_alignment_base_dir, str(i))
            for i in xrange(sys.maxint)
            if not os.path.exists(
                    os.path.join(contig_alignment_base_dir, str(i)))).next()

    os.mkdir(contig_alignment_dir)

    # Concatenate contig fastas for alignment
    contig_fastas = [get_fasta(c) for c in contig_list]
    contig_concat = os.path.join(contig_alignment_dir, 'contig_concat.fa')
    with open(contig_concat, 'w') as fh:
        subprocess.check_call(' '.join(['cat'] + contig_fastas),
                shell=True, executable=settings.BASH_PATH, stdout=fh)

    # Write mobile elements to file
    ref_genome = sample_alignment.alignment_group.reference_genome
    genbank_query = ref_genome.dataset_set.filter(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK)

    use_me_alignment = bool(genbank_query.count())
    if use_me_alignment:
        genbank_path = genbank_query[0].get_absolute_location()
        me_concat_fasta = os.path.join(contig_alignment_dir, 'me_concat.fa')
        if not sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.MOBILE_ELEMENT_FASTA):
            write_me_features_multifasta(genbank_path, me_concat_fasta)
            add_dataset_to_entity(
                    sample_alignment,
                    Dataset.TYPE.MOBILE_ELEMENT_FASTA,
                    Dataset.TYPE.MOBILE_ELEMENT_FASTA,
                    me_concat_fasta)

        contig_alignment_to_me_bam = os.path.join(
                contig_alignment_dir, 'contig_alignment_to_me.bam')

        me_concat_fasta = sample_alignment.dataset_set.get(
                type=Dataset.TYPE.MOBILE_ELEMENT_FASTA).get_absolute_location()

        if not os.path.exists(contig_alignment_to_me_bam):
            ensure_bwa_index(me_concat_fasta)
            simple_align_with_bwa_mem(
                    contig_concat,
                    me_concat_fasta,
                    contig_alignment_to_me_bam)

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
    if use_me_alignment:
        add_me_alignment_to_graph(G, contig_alignment_to_me_bam)


    # Add SEQUENCE_GRAPH_PICKLE dataset to sample alignment
    graph_pickle_path = os.path.join(
            contig_alignment_dir,
            'sequence_graph.pickle')
    nx.write_gpickle(G, graph_pickle_path)
    add_dataset_to_entity(
            sample_alignment,
            Dataset.TYPE.SEQUENCE_GRAPH_PICKLE,
            Dataset.TYPE.SEQUENCE_GRAPH_PICKLE,
            graph_pickle_path)

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
        var_dict_list = [parse_path_into_ref_alt(iv_pair, contig_qname_to_uid,
                sample_alignment)
                for iv_pair in trans_iv_pairs]

        var_dict_list = [var_d for var_d in var_dict_list
                if any([var_d['ref_seq'], var_d['alt_seq']])]

        if use_me_alignment:
            me_trans_iv_pairs = me_translocation_walk(G)

            me_var_dict_list = [parse_path_into_ref_alt(
                    iv_pair, contig_qname_to_uid,
                    sample_alignment)
                for iv_pair in me_trans_iv_pairs]
        else:
            me_var_dict_list = []

    else:
        print 'Translocation walk not implemented for multi-chromosomal refs'
        var_dict_list = []
        me_var_dict_list = []

    return placeable_contigs, var_dict_list, me_var_dict_list


def avg_coverage(bam, chromosome, start, end):
    """Average read coverage between start and end positions on
    the chromosome as indicated by the bam file whose path is passed as an arg
    """

    alignment_file = pysam.AlignmentFile(bam)
    coverage_atcg = alignment_file.count_coverage(
            str(chromosome), start, end)

    avg_coverage = sum(max(atcg) for atcg in zip(*[coverage_atcg[i] for i in
            range(4)])) / (end - start)

    return avg_coverage


def add_me_alignment_to_graph(G, contig_alignment_bam, add_rc_me_seqs=True):
    """Add (sequence vertex, sequence vertex) edges between contigs and
    mobile elements to DiGraph G as indicated by the contig_alignment_bam
    alignment file
    """

    contigs_intervals = G.contig_intervals_list
    me_interval_dict = {}

    contig_alignmentfile = pysam.AlignmentFile(contig_alignment_bam)
    for me_ref_name, me_len in zip(
            contig_alignmentfile.references, contig_alignmentfile.lengths):
        # Create sequence interval instances for each ME
        me_interval_dict[me_ref_name] = SequenceIntervals(
                me_ref_name, me_len)

    # Iterate over aligned contig 'reads' in contig alignment to ref bam
    for read in contig_alignmentfile:

        contig_intervals = contigs_intervals[read.qname]

        match_regions = get_match_regions(read)
        for match_region in match_regions:
            if match_region.read_start != 0:
                # Insert break into sequence intervals for ref and contig
                ref_name = contig_alignmentfile.getrname(read.reference_id)
                me_seq_intervals = me_interval_dict[ref_name]
                me_vert = me_seq_intervals.insert_vertex(
                        match_region.ref_start)
                contig_vert = contig_intervals.insert_vertex(
                        match_region.read_start)

                # Add an edge to G for the junction
                G.add_edge(contig_vert, me_vert,
                        weight=match_region.length, is_rc=read.is_reverse,
                        match_region=match_region)
                G.add_edge(me_vert, contig_vert,
                        weight=match_region.length, is_rc=read.is_reverse,
                        match_region=match_region)

            if match_region.read_end != contig_intervals.length:
                # Insert break into sequence intervals for ref and contig
                ref_name = contig_alignmentfile.getrname(read.reference_id)
                me_seq_intervals = me_interval_dict[ref_name]
                me_vert = me_seq_intervals.insert_vertex(
                        match_region.ref_end)
                contig_vert = contig_intervals.insert_vertex(
                        match_region.read_end)

                # Add an edge to G for the junction
                G.add_edge(me_vert, contig_vert,
                        weight=match_region.length, is_rc=read.is_reverse,
                        match_region=match_region)
                G.add_edge(contig_vert, me_vert,
                        weight=match_region.length, is_rc=read.is_reverse,
                        match_region=match_region)

    G.me_interval_dict = me_interval_dict

    if add_rc_me_seqs:
        add_rc_of_me_to_graph(G)

    # Add inter-contig sequence edges
    for contig_intervals in contigs_intervals.values():
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


def add_rc_of_me_to_graph(G):
    """Creates a SequenceIntervals instance for the reverse complement of each
    mobile element and adds the appropriate edges to DiGraph G
    """

    me_interval_dict = G.me_interval_dict

    rc_interval_dict = {}
    for me_interval in me_interval_dict.values():

        l = me_interval.length
        seq_uid = me_interval.seq_uid + '_RC'
        rc_interval = SequenceIntervals(seq_uid, l)

        for v in me_interval.vertices:
            new_v = rc_interval.insert_vertex(l - v.pos)
            for _, endpoint, data in G.edges([v], data=True):
                contig_interval = endpoint.parent
                contig_l = contig_interval.length
                new_contig_v = contig_interval.insert_vertex(
                        contig_l - endpoint.pos)
                data['is_rc'] = not data['is_rc']
                G.add_edge(new_v, new_contig_v, **data)
                G.add_edge(new_contig_v, new_v, **data)

        rc_interval_dict[seq_uid] = rc_interval

    G.me_interval_dict.update(rc_interval_dict)


def add_alignment_to_graph(G, contig_alignment_bam):
    """Add (sequence vertex, sequence vertex) edges between contigs and
    the reference to DiGraph G as indicated by the contig_alignment_bam
    alignment file
    """

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

            if match_region.read_start != 0:
                # Insert break into sequence intervals for ref and contig
                ref_vert = ref_intervals.insert_vertex(match_region.ref_start)
                contig_vert = contig_intervals.insert_vertex(
                        match_region.read_start)

                # Add an edge to G for the junction
                G.add_edge(contig_vert, ref_vert,
                        weight=match_region.length, is_rc=read.is_reverse,
                        match_region=match_region)

            if match_region.read_end != contig_intervals.length:
                # Insert break into sequence intervals for ref and contig
                ref_vert = ref_intervals.insert_vertex(match_region.ref_end)
                contig_vert = contig_intervals.insert_vertex(
                        match_region.read_end)

                # Add an edge to G for the junction
                G.add_edge(ref_vert, contig_vert,
                        weight=match_region.length, is_rc=read.is_reverse,
                        match_region=match_region)

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
    """Walk the graph and return InsertionVertices objects corresponding to
    ref - contig - contig - ref paths that represent novel sequence insertions
    """

    ref_seq_uid = G.ref_intervals.vertices[0].seq_uid
    contig_seq_uid_set = set(ci.seq_uid for ci in
            G.contig_intervals_list.values())

    def ref_neighbors(vert):
        return [v for v in G.neighbors(vert) if v.seq_uid == ref_seq_uid]

    def contig_neighbors(vert):
        return [v for v in G.neighbors(vert) if
                v.seq_uid in contig_seq_uid_set]

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

                extend_items = []
                for n in contig_neighbors(exit_contig):
                    if n not in visited:
                        extend_items.append(n)

                queue.extend(extend_items)

    return iv_list


def me_translocation_walk(G):
    """Walk the graph and return InsertionVertices tuples corresponding to
    ref - contig - contig - ME - ME - contig - contig - ref paths that
    represent mobile element insertions

    Args:
        G: networkx.DiGraph representation of sequence graph

    Returns:
        filtered: (InsertionVertices, InsertionVertices) tuples corresponding
            to ME insertion paths through the sequence graph
    """

    assert hasattr(G, 'me_interval_dict')

    ref_seq_uid_set = set([G.ref_intervals.seq_uid])

    me_seq_uid_set = set(si.seq_uid for si in G.me_interval_dict.values())
    ref_seq_uid_set = ref_seq_uid_set | me_seq_uid_set

    contig_seq_uid_set = set(ci.seq_uid for ci in
            G.contig_intervals_list.values())

    def ref_neighbors(vert):
        if vert not in G:
            return []
        return [v for v in G.neighbors(vert) if v.seq_uid in ref_seq_uid_set]

    def contig_neighbors(vert):
        if vert not in G:
            return []
        return [v for v in G.neighbors(vert) if v.seq_uid in
                contig_seq_uid_set]

    forward_edges = []
    back_edges = []

    me_vertices = reduce(lambda x, y: x + y,
            [si.vertices for si in G.me_interval_dict.values()])

    dset = set()
    for exit_ref in G.ref_intervals.vertices + me_vertices:
        for enter_contig in contig_neighbors(exit_ref):
            queue = set([enter_contig])
            visited = []
            while queue:
                exit_contig = queue.pop()
                for enter_ref in ref_neighbors(exit_contig):

                    iv = InsertionVertices(
                            exit_ref, enter_contig, exit_contig,
                            enter_ref)
                    dset.add(iv)

                    if exit_ref.pos < enter_ref.pos:
                        forward_edges.append(iv)
                    else:
                        back_edges.append(iv)

                visited.append(exit_contig)
                queue.update([n for n in contig_neighbors(exit_contig)
                        if n not in visited])

    sorted_by_exit_ref = sorted(forward_edges + back_edges,
            key=lambda x: x.exit_ref.pos)

    sorted_by_enter_ref = sorted(forward_edges + back_edges,
            key=lambda x: x.enter_ref.pos)

    me_iv_pairs = []
    i = 0

    me_sorted_by_exit_ref = [e for e in sorted_by_exit_ref if (
            e.enter_ref.seq_uid in me_seq_uid_set and
            e.exit_ref.seq_uid == G.ref_intervals.seq_uid)]

    me_sorted_by_enter_ref = [e for e in sorted_by_enter_ref if (
            e.exit_ref.seq_uid in me_seq_uid_set and
            e.enter_ref.seq_uid == G.ref_intervals.seq_uid)]

    OVERLAP_TOLERANCE = 400
    for enter_iv in me_sorted_by_exit_ref:

        while (me_sorted_by_enter_ref[i].enter_ref.pos <
                enter_iv.exit_ref.pos - OVERLAP_TOLERANCE and
                i < len(me_sorted_by_enter_ref) - 1):
            i += 1

        j = i
        exit_iv = me_sorted_by_enter_ref[j]
        deletion = exit_iv.enter_ref.pos - enter_iv.exit_ref.pos
        while deletion < MAX_TRANS_DELETION:

            # Add to list if path enters and exits from the same ME
            if enter_iv.enter_ref.seq_uid == exit_iv.exit_ref.seq_uid:
                me_iv_pairs.append((enter_iv, exit_iv))

            if j == len(me_sorted_by_enter_ref) - 1:
                break

            j += 1
            exit_iv = me_sorted_by_enter_ref[j]
            deletion = exit_iv.enter_ref.pos - enter_iv.exit_ref.pos

    filtered = strand_filter(G, me_iv_pairs)
    filtered = me_length_filter(filtered)
    filtered = match_region_filter(G, filtered)

    return filtered


def strand_filter(G, me_iv_pairs):
    """Returns list of (InsertionVertices, InsertionVertices) tuples
    with those containing paths of inconsistent polarity filtered out
    """

    def _is_edge_rc(u, v):
        data = G.get_edge_data(u, v)
        return data['is_rc']

    filtered = []
    for iv_pair in me_iv_pairs:

        enter_iv, exit_iv = iv_pair

        r_to_c_strand = _is_edge_rc(enter_iv.exit_ref, enter_iv.enter_contig)
        c_to_m_strand = _is_edge_rc(enter_iv.exit_contig, enter_iv.enter_ref)
        m_to_c_strand = _is_edge_rc(exit_iv.exit_ref, exit_iv.enter_contig)
        c_to_r_strand = _is_edge_rc(exit_iv.exit_contig, exit_iv.enter_ref)

        enter_me_agreement = r_to_c_strand == c_to_m_strand
        exit_me_agreement = m_to_c_strand == c_to_r_strand

        me_strand_agreement = (
                enter_iv.enter_ref.seq_uid.endswith('_RC') ==
                exit_iv.exit_ref.seq_uid.endswith('_RC'))

        if me_strand_agreement and enter_me_agreement and exit_me_agreement:
            filtered.append(iv_pair)

    return filtered


def me_length_filter(me_iv_pairs, min_length=100):
    """Returns list of (InsertionVertices, InsertionVertices) tuples
    with those containing paths going backwards through the ME sequence
    filtered out
    """

    filtered = []
    for iv_pair in me_iv_pairs:
        enter_iv, exit_iv = iv_pair
        me_seq_len = exit_iv.exit_ref.pos - enter_iv.enter_ref.pos

        if me_seq_len > min_length:
            filtered.append(iv_pair)

    return filtered


def match_region_filter(G, me_iv_pairs):
    """Returns list of (InsertionVertices, InsertionVertices) tuples,
    filtering out paths corresponding to null variants and only keeping
    a single copy of synonymous variants
    """

    def _me_seq_length(iv_pair):
        enter_iv, exit_iv = iv_pair
        return exit_iv.exit_ref.pos - enter_iv.enter_ref.pos

    filtered = {}
    for iv_pair in me_iv_pairs:

        enter_iv, exit_iv = iv_pair

        enter_me_edge = (enter_iv.exit_contig, enter_iv.enter_ref)
        exit_me_edge = (exit_iv.exit_ref, exit_iv.enter_contig)

        enter_me_match_region = G.get_edge_data(*enter_me_edge)['match_region']
        exit_me_match_region = G.get_edge_data(*exit_me_edge)['match_region']

        if enter_me_match_region == exit_me_match_region:
            continue

        key = (enter_iv.exit_ref, exit_iv.enter_ref,
                enter_me_match_region, exit_me_match_region)

        redundant = filtered.get(key, None)
        if not redundant:
            filtered[key] = iv_pair
        else:
            if _me_seq_length(iv_pair) > _me_seq_length(redundant):
                filtered[key] = iv_pair

    return filtered.values()


def translocation_walk(G):
    """Walk the graph and return InsertionVertices tuples corresponding to
    ref - contig - contig - ref - ref - contig - contig - ref paths that
    represent translocations within the reference

    Args:
        G: networkx.DiGraph representation of sequence graph

    Returns:
        filtered: (InsertionVertices, InsertionVertices) tuples corresponding
            to translocation paths through the sequence graph
    """

    ref_seq_uid_set = set([G.ref_intervals.seq_uid])

    contig_seq_uid_set = set(ci.seq_uid for ci in
            G.contig_intervals_list.values())

    def ref_neighbors(vert):
        if vert not in G:
            return []
        return [v for v in G.neighbors(vert) if v.seq_uid in ref_seq_uid_set]

    def contig_neighbors(vert):
        if vert not in G:
            return []
        return [v for v in G.neighbors(vert) if v.seq_uid in
                contig_seq_uid_set]

    forward_edges = []
    back_edges = []


    dset = set()
    for exit_ref in G.ref_intervals.vertices:
        for enter_contig in contig_neighbors(exit_ref):
            queue = set([enter_contig])
            visited = []
            while queue:
                exit_contig = queue.pop()
                for enter_ref in ref_neighbors(exit_contig):

                    iv = InsertionVertices(
                            exit_ref, enter_contig, exit_contig,
                            enter_ref)
                    dset.add(iv)

                    if exit_ref.pos < enter_ref.pos:
                        forward_edges.append(iv)
                    else:
                        back_edges.append(iv)

                visited.append(exit_contig)
                queue.update([n for n in contig_neighbors(exit_contig)
                        if n not in visited])

    sorted_by_exit_ref = sorted(forward_edges + back_edges,
            key=lambda x: x.exit_ref.pos)

    sorted_by_enter_ref = sorted(forward_edges + back_edges,
            key=lambda x: x.enter_ref.pos)

    iv_pairs = []
    i = 0

    OVERLAP_TOLERANCE = 400
    for enter_iv in sorted_by_exit_ref:

        while (sorted_by_enter_ref[i].enter_ref.pos <
                enter_iv.exit_ref.pos - OVERLAP_TOLERANCE and
                i < len(sorted_by_enter_ref) - 1):
            i += 1

        j = i
        exit_iv = sorted_by_enter_ref[j]
        deletion = exit_iv.enter_ref.pos - enter_iv.exit_ref.pos
        while deletion < MAX_TRANS_DELETION:

            # Length of translocation sequence
            trans_length = exit_iv.exit_ref.pos - enter_iv.enter_ref.pos

            if (all(v.seq_uid in ref_seq_uid_set for v in
                    [exit_iv.exit_ref, exit_iv.enter_ref,
                    enter_iv.exit_ref, enter_iv.enter_ref]) and (
                    MIN_TRANS_LENGTH < trans_length < MAX_TRANS_LENGTH)):
                iv_pairs.append((enter_iv, exit_iv))

            if j == len(sorted_by_enter_ref) - 1:
                break

            j += 1
            exit_iv = sorted_by_enter_ref[j]
            deletion = exit_iv.enter_ref.pos - enter_iv.exit_ref.pos

    filtered = match_region_filter(G, iv_pairs)

    return filtered


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
            raise TypeError
        if rhs.seq_uid != self.seq_uid:
            raise ValueError

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

    def blank_copy(self):
        return SequenceIntervals(self.seq_uid, self.length, self.tag)


def get_fasta(has_fasta):
    return has_fasta.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()


MatchRegion = namedtuple('MatchRegion',
            ['ref_start', 'ref_end', 'read_start', 'read_end', 'length'])


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
    match_regions = []
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


def parse_path_into_ref_alt(path_list, contig_qname_to_uid,
        sample_alignment):
    """Takes the (InsertionVertices, Insertion Vertices) tuple path_list
    and returns a variant in the form of a dictionary with keys:
    chromosome, pos, ref_seq, alt_seq
    """

    ref_genome = sample_alignment.alignment_group.reference_genome
    ref_uid = ref_genome.uid
    ref_fasta = get_fasta(ref_genome)
    with open(ref_fasta) as fh:
        ref_seqrecord = SeqIO.parse(fh, 'fasta').next()
        ref_seq = str(ref_seqrecord.seq)
        ref_chromosome = ref_seqrecord.id

    def _seq_str(enter_vert, exit_vert):
        if enter_vert.seq_uid == ref_uid:
            return ref_seq[enter_vert.pos: exit_vert.pos]

        if enter_vert.seq_uid.startswith('ME_'):

            rc = enter_vert.seq_uid[-3:] == '_RC'
            if rc:
                seq_uid = enter_vert.seq_uid[:-3]
            else:
                seq_uid = enter_vert.seq_uid

            me_fasta = sample_alignment.dataset_set.get(
                    type=Dataset.TYPE.MOBILE_ELEMENT_FASTA
            ).get_absolute_location()

            seq_rec = (sr for sr in SeqIO.parse(me_fasta, 'fasta')
                    if sr.id == seq_uid).next()

            assert enter_vert.pos <= exit_vert.pos
            if rc:
                return str(seq_rec.seq.reverse_complement()[
                        enter_vert.pos: exit_vert.pos])
            else:
                return str(seq_rec.seq[enter_vert.pos: exit_vert.pos])

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
            yield path_list_concat[i:i + 2]
            i += 2

    seq_list = []
    for enter_vert, exit_vert in _seq_interval_iter():

        seq_len = exit_vert.pos - enter_vert.pos
        seq_uid = enter_vert.seq_uid
        if seq_uid.startswith('ME_'):
            if seq_len == 0:
                seq_list.append('')
            else:
                seq_list.append(_seq_str(enter_vert, exit_vert))

        elif seq_len < 0:
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

    var_dict = {
        'chromosome': ref_chromosome,
        'pos': ref_start,
        'ref_seq': ref_seq,
        'alt_seq': alt_seq
    }

    translocation_seq_verts = [path_list_concat[3], path_list_concat[4]]
    seq_names = []
    me_prefix = 'ME_insertion_sequence:'
    is_me = False
    for v in translocation_seq_verts:
        name = v.seq_uid
        if name.startswith(me_prefix):
            is_me = True
            name = name[len(me_prefix):]
        seq_names.append(name)

    if is_me:
        # Assert the two vertices of the translocation in the path
        # are on the same mobile element sequence
        assert len(set(seq_names)) == 1

        name = seq_names[0]
        is_rc = name[-3:] == '_RC'
        if is_rc:
            name = name[:-3]

        var_dict['MEINFO'] = {
                'name': seq_names[0],
                'start': translocation_seq_verts[0].pos,
                'end': translocation_seq_verts[1].pos,
                'polarity': '-' if is_rc else '+'
        }

    return var_dict
