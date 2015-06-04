#!/usr/bin/env python

from optparse import OptionParser
import sys

from sam_parse import extractCigarOps
from sam_parse import SAM


def get_clipped_reads(inFile, outFile=None, clipping_threshold=8):

    if outFile:
        out_data = open(outFile, 'w+')

    def write_out(line):
        if outFile:
            out_data.write(line)
        else:
            print line

    if inFile == "stdin":
        data = sys.stdin
    else:
        data = open(inFile, 'r')
    for line in data:
        split = 0
        if line[0] == '@':
            print line.strip()
            continue
        samList = line.strip().split('\t')
        sam = SAM(samList)
        for el in sam.tags:
            if "SA:" in el:
                split = 1

        if not split:  # So as to not repeat reads from extractSplitReads
            samList[0] = sam.query
            readCigar = sam.cigar
            readCigarOps = extractCigarOps(readCigar, sam.flag)

            # readCigarOps is a list of the form [cigar, cigar, cigar]
            # where cigar.op and cigar.length are the letter and number
            for cigar_op in readCigarOps:
                if cigar_op.op in ['H', 'S'] and (
                            cigar_op.length >= clipping_threshold):
                    print "\t".join(samList)
                    continue


def main():

    usage = """%prog -i <file>"""

    parser = OptionParser(usage)
    parser.add_option(
        "-i", "--inFile", dest="inFile",
        help="A SAM file or standard input (-i stdin).",
        metavar="FILE")
    (opts, args) = parser.parse_args()
    if opts.inFile is None:
        parser.print_help()
        print
    else:
        try:
            get_clipped_reads(opts.inFile)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return
if __name__ == "__main__":
    sys.exit(main())
