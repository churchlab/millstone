#!/usr/bin/env python

from optparse import OptionParser
import pickle
import sys

from sam_parse import extractCigarOps
from sam_parse import SAM


def count_matches(inFile):
    """
    Input: SAM file
    Output: Dictionary with QNAME keys and lists of
            lengths of consecutive base-pair matching runs as values
    """
    if inFile == "stdin":
        data = sys.stdin
    else:
        data = open(inFile, 'r')
    match_counts = {}
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
        match_count = 0
        for cigar_op in readCigarOps:
            if cigar_op.op == 'M':
                match_count += cigar_op.length

        if sam.query in match_counts:
            match_counts[sam.query].append(match_count)
        else:
            match_counts[sam.query] = [match_count]

    return match_counts


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
            count_dict = count_matches(opts.inFile)
            print pickle.dumps(count_dict)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return
if __name__ == "__main__":
    sys.exit(main())
