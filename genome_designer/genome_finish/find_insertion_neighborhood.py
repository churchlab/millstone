#!/usr/bin/env python

from optparse import OptionParser
import pickle
import sys

from sam_parse import calcQueryPosFromCigar
from sam_parse import extractCigarOps
from sam_parse import SAM


def find_insertion_neighborhood(inFile):
    """
    Input: SAM file
    Output: Dictionary with QNAME keys and lists of
            lengths of consecutive base-pair matching runs as values
    """
    if inFile == "stdin":
        data = sys.stdin
    else:
        data = open(inFile, 'r')
    left_matches = 0
    right_matches = 0
    # contig_length = 0
    left_end = -1
    right_end = -1
    ref = None
    for line in data:
        if line[0] == '@':
            continue
        samList = line.strip().split('\t')
        sam = SAM(samList)
        samList[0] = sam.query
        readCigar = sam.cigar
        readCigarOps = extractCigarOps(readCigar, sam.flag)

        # readCigarOps is a list of the form [cigar, cigar, cigar]
        # where cigar.op and cigar.length are the letter and number
        current_match_stretch = 0
        following_clipping = False
        for i, cigar_op in enumerate(readCigarOps):
            if cigar_op.op == 'M':
                current_match_stretch += cigar_op.length
            elif cigar_op.op in ['H', 'S']:
                if current_match_stretch > left_matches:
                    left_matches = current_match_stretch
                    left_end = sam.pos
                    chromosome = sam.ref
                if (following_clipping and
                        current_match_stretch > right_matches):
                    right_matches = current_match_stretch
                    queryPos = calcQueryPosFromCigar(readCigarOps)
                    alignment_length = queryPos.qePos - queryPos.qsPos
                    # assert 1==0, 'qsPos:' + str(queryPos.qsPos) + ' qePos:' + str(queryPos.qePos) + ' qLen:' + str(queryPos.qLen)
                    right_end = sam.pos + alignment_length

                following_clipping = True
                current_match_stretch = 0

        if following_clipping and current_match_stretch > right_matches:
            right_matches = current_match_stretch
            queryPos = calcQueryPosFromCigar(readCigarOps)
            alignment_length = queryPos.qePos - queryPos.qsPos
            # assert 1==0, 'sam.pos:' + str(sam.pos) + ' qsPos:' + str(queryPos.qsPos) + ' qePos:' + str(queryPos.qePos) + ' qLen:' + str(queryPos.qLen)
            right_end = sam.pos + alignment_length

    assert left_end != -1 and right_end != -1, ('left or right end not found;'
            'left end:', left_end, 'right end:', right_end)

    assert left_end < right_end, ('right end comes before left end;',
            'left end:', left_end, 'right end:', right_end)

    locations_dict = {
        'chromosome_seqrecord_id': chromosome,
        'left_end': left_end,
        'left_matches': left_matches,
        'right_end': right_end,
        'right_matches': right_matches
    }
    return locations_dict


def main():

    usage = """%prog -i <file>"""

    parser = OptionParser(usage)
    parser.add_option(
            "-i", "--inFile", dest="inFile",
            help="A SAM file or standard input (-i stdin).", metavar="FILE")
    (opts, args) = parser.parse_args()
    if opts.inFile is None:
        parser.print_help()
        print
    else:
        try:
            locations_dict = find_insertion_neighborhood(opts.inFile)
            print pickle.dumps(locations_dict)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return
if __name__ == "__main__":
    sys.exit(main())
