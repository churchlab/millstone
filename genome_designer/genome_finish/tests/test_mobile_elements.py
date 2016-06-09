import networkx as nx
import os

from django.conf import settings
from django.test import TestCase

from genome_finish.graph_contig_placement import add_alignment_to_graph
from genome_finish.graph_contig_placement import add_me_alignment_to_graph
from genome_finish.graph_contig_placement import me_translocation_walk
from genome_finish.graph_contig_placement import SequenceIntervals

from main.testing_util import create_common_entities

GF_TEST_DIR = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test')


class TestMobileElements(TestCase):
    """
    This is Line 5 from Tenaillon 2013.
    This test should give two IS150 elements at 1305489 and 968327
    and an IS186 element at 4517605 (which is in the Ancestor also.)
    """

    def setUp(self):
        self.common_data = create_common_entities()
        self.project = self.common_data['project']

        self.pos_elem_tups = [
                    (968327, 'IS150'),
                    (1305482, 'IS150'),
                    (4517605, 'IS186')]

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

        self.me_iv_pairs = me_translocation_walk(G)


    def test_3_mobile_elements(self):
        """
        Assert elements are found in the correct locations.
        """

        def _is_correct(iv_pair, pos, element_name='', tolerance=100):
            if abs(iv_pair[0][0].pos - pos) > tolerance:
                return False

            found_element_name = iv_pair[0][3].seq_uid
            me_prefix = 'ME_insertion_sequence:'
            if not found_element_name.startswith(me_prefix + element_name):
                return False

        return True

        # Assert all expected MEs are found
        for pos, elem in self.pos_elem_tups:
            self.assertTrue(any(_is_correct(iv_pair, pos, elem)
                    for iv_pair in self.me_iv_pairs))

        # Assert no redundant MEs
        for pos, elem in self.pos_elem_tups:
            self.assertTrue(sum(_is_correct(iv_pair, pos, elem)
                    for iv_pair in self.me_iv_pairs) == 1)

        # Assert no false positive MEs
        assert len(self.me_iv_pairs) == len(self.pos_elem_tups)
