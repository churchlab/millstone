#!/usr/bin/env python

from optparse import OptionParser
import pickle
import sys

from sam_parse import calcQueryPosFromCigar
from sam_parse import cigarLength
from sam_parse import extractCigarOps
from sam_parse import SAM

INTER_CLIPPED_END_DISTANCE_TOLERANCE = 8


def find_insertion_neighborhood(inFile):
    """
    Input: SAM file
    Output: Dictionary with insertion region positional information
    """
    if inFile == "stdin":
        data = sys.stdin
    else:
        data = open(inFile, 'r')

    left_sam = None
    right_sam = None

    for line in data:
        # Skip header
        if line[0] == '@':
            continue
        samList = line.strip().split('\t')
        sam = SAM(samList)
        readCigarOps = extractCigarOps(sam.cigar, sam.flag)

        # Grab right and left alignments by looking the matching end
        if readCigarOps[0].op not in ['H', 'S']:
            left_sam = sam
            left_readCigarOps = readCigarOps
        elif readCigarOps[-1] not in ['H', 'S']:
            right_sam = sam
            right_readCigarOps = readCigarOps

    if left_sam is None:
        error_string = 'Expected homology on left end of contig not found'
    elif right_sam is None:
        error_string = 'Expected homology on right end of contig not found'
    else:
        # Find left and right alignment regions
        left_cigar_ops = [s.op for s in left_readCigarOps]
        right_cigar_ops = [s.op for s in right_readCigarOps]

        if 'H' not in left_cigar_ops and 'S' not in left_cigar_ops:
            error_string = ('Left side of split alignment not clipped, ' +
                'indicating a lack of evidence for the insertion')
        elif 'H' not in right_cigar_ops and 'S' not in right_cigar_ops:
            error_string = ('Right side of split alignment not clipped, ' +
                'indicating a lack of evidence for the insertion')
        else:
            left_start_pos = left_sam.pos
            left_clip_index = min(
                    left_cigar_ops.index('H')
                    if 'H' in left_cigar_ops else float('inf'),
                    left_cigar_ops.index('S')
                    if 'S' in left_cigar_ops else float('inf'))
            left_end_pos = left_sam.pos + (
                    cigarLength(left_readCigarOps[:left_clip_index]))

            right_clip_index = (len(right_cigar_ops) - 1 - min(
                    list(reversed(right_cigar_ops)).index('H')
                    if 'H' in right_cigar_ops else float('inf'),
                    list(reversed(right_cigar_ops)).index('S')
                    if 'S' in right_cigar_ops else float('inf')))
            right_start_pos = right_sam.pos + (
                    cigarLength(right_readCigarOps[:right_clip_index+1]))
            right_end_pos = right_sam.pos + (cigarLength(right_readCigarOps))

            # Find any potential errors with posited insertion region
            if left_sam.ref != right_sam.ref:
                error_string = ('Left and right ends of alignment map to ' +
                    'different chromosomes')
            if left_start_pos > right_start_pos:
                error_string = ('Left half of split alignment is clipped on ' +
                        'left end rather than right, indicating a more ' +
                        'complex structural variant')
            elif (abs(left_end_pos - right_start_pos) >
                    INTER_CLIPPED_END_DISTANCE_TOLERANCE):
                error_string = ('Distance between left and right clipped ' +
                        'ends of the split alignment too far apart, ' +
                        'indicating a repeat or deletion in the reference')
            elif (right_end_pos - left_start_pos >
                    calcQueryPosFromCigar(left_readCigarOps).qLen):
                error_string = ('Contig shorter than distance between ' +
                        'beginning and end of split alignment, perhaps ' +
                        'indicating a deletion in the reference')
            else:
                locations_dict = {
                    'chromosome_seqrecord_id': left_sam.ref,
                    'left_end': left_start_pos,
                    'right_end': right_end_pos,
                }
                return locations_dict

    return {'error_string': error_string}


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
