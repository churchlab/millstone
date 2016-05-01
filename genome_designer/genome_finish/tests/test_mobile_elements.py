import networkx as nx
import os

from django.conf import settings
from django.test import TestCase

from genome_finish.graph_contig_placement import add_alignment_to_graph
from genome_finish.graph_contig_placement import add_me_alignment_to_graph
from genome_finish.graph_contig_placement import me_translocation_walk
from genome_finish.graph_contig_placement import SequenceIntervals

GF_TEST_DIR = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test')


class TestMobileElements(TestCase):

    def test_2_mobile_elements(self):

        graph_test_dir = os.path.join(GF_TEST_DIR, 'graph_tests')

        contig_alignment_bam = os.path.join(
                graph_test_dir, 'contig_alignment.bam')

        contig_alignment_to_me_bam = os.path.join(
                graph_test_dir, 'contig_alignment_to_me.bam')

        # Create graph
        G = nx.DiGraph()

        # Create sequence interval instances for reference and each contig
        ref_intervals = SequenceIntervals(
                'ref', 4622887, tag='ref')

        G.ref_intervals = ref_intervals

        add_alignment_to_graph(G, contig_alignment_bam)
        add_me_alignment_to_graph(G, contig_alignment_to_me_bam)

        me_iv_pairs = me_translocation_walk(G)

        def _is_correct(iv_pair, pos, element_name='', tolerance=100):
            if abs(iv_pair[0][0].pos - pos) > tolerance:
                return False

            found_element_name = iv_pair[0][3].seq_uid
            me_prefix = 'ME_insertion_sequence:'
            if not found_element_name.startswith(me_prefix + element_name):
                return False

            return True

        pos_elem_tups = [(968327, 'IS150'), (1305482, 'IS150')]

        # Assert all expected ME's are found
        for pos, elem in pos_elem_tups:
            self.assertTrue(any(_is_correct(iv_pair, pos, elem)
                    for iv_pair in me_iv_pairs))

        # Assert no redundant ME's
        for pos, elem in pos_elem_tups:
            self.assertTrue(sum(_is_correct(iv_pair, pos, elem)
                    for iv_pair in me_iv_pairs) == 1)

        # Assert no false positive ME's
        assert len(me_iv_pairs) == len(pos_elem_tups)
