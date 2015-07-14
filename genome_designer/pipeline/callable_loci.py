import pysam
import itertools
import sys
from collections import defaultdict

from django.conf import settings

MIN_MAPQ = settings.CL__MIN_MAPQ
MAX_DEPTH = settings.CL__MAX_DEPTH
MIN_DEPTH = settings.CL__MIN_DEPTH
MIN_LOWMAPQ_DEPTH = settings.CL__MIN_LOWMAPQ_DEPTH
MAX_LOWMAP_FRAC = settings.CL__MAX_LOWMAP_FRAC
MERGE_DIST = settings.CL__MERGE_DIST

def get_callable_loci(
        bam_filename,
        bed_output,
        chrom=None,
        start=None,
        end=None):

    # all bed lines, grouped by flag
    bed_lines = defaultdict(list)

    def _save_bed_line(chrom, bed_start, bed_end, flag):
        if not flag: return

        new_bed_line = {
            'chrom': chrom,
            'start': bed_start,
            'end': bed_end,
            'name': flag
        }

        # merge with previous if within MERGE_DIST:
        try:
            prev_bed_line = bed_lines[flag][-1]

            same_chrom = prev_bed_line['chrom'] == chrom
            within_merge = prev_bed_line['end'] + MERGE_DIST > bed_start

            if same_chrom and within_merge:
                prev_bed_line['end'] = bed_end
            else:
                bed_lines[flag].append(new_bed_line)

        # can't merge if this is the first bed line
        except IndexError:
            bed_lines[flag].append(new_bed_line)

    bamfile = pysam.AlignmentFile(bam_filename, 'rb')

    chrom_list = bamfile.references
    chrom_lens = bamfile.lengths

    if chrom:
        chrom_dict = dict(zip(chrom_list, chrom_lens))
        chrom_list = [chrom]
        if not start:
            c_starts = [0]
        else:
            c_starts = [start]
        if not end:
            c_ends = chrom_dict[chrom]
        else:
            c_ends = [end]
    else:
        c_starts = [0]*len(chrom_list)
        c_ends = chrom_lens

    for chrom, c_start, c_end in zip(chrom_list, c_starts, c_ends):
        curr_flag = None
        next_pos = c_start
        bed_start = c_start

        for pileup_col in bamfile.pileup(chrom,
                start=c_start, end=c_end, truncate=True):

            depth = pileup_col.nsegments
            badmapq = 0
            altaligns = 0

            for p in pileup_col.pileups:
                a = p.alignment
                tag_as = a.get_tag('AS')
                tag_xs = a.get_tag('XS')
                badmapq += a.mapping_quality < MIN_MAPQ
                altaligns += tag_as <= tag_xs

            # NO_COVERAGE
            if next_pos != pileup_col.pos:
                if curr_flag != 'NO_COVERAGE':
                    _save_bed_line(chrom, bed_start, next_pos-1, curr_flag)
                _save_bed_line(chrom, next_pos, pileup_col.pos-1, 'NO_COVERAGE')
                curr_flag = None

            # LOW_COVERAGE
            if depth < MIN_DEPTH:
                if curr_flag != 'LOW_COVERAGE':
                    _save_bed_line(
                            chrom, bed_start, pileup_col.pos-1, curr_flag)
                    bed_start = pileup_col.pos
                    curr_flag = 'LOW_COVERAGE'

            if depth >= MIN_LOWMAPQ_DEPTH:
                # POOR_MAPQ
                if badmapq >= MAX_LOWMAP_FRAC*depth:
                    if curr_flag != 'POOR_MAP_QUALITY':
                        _save_bed_line(
                                chrom, bed_start, pileup_col.pos-1, curr_flag)
                        bed_start = pileup_col.pos
                        curr_flag = 'POOR_MAP_QUALITY'

                # NONUNIQUE
                elif altaligns >= MAX_LOWMAP_FRAC*depth:
                    if curr_flag != 'NONUNIQUE_ALIGNMENTS':
                        _save_bed_line(
                                chrom, bed_start, pileup_col.pos-1, curr_flag)
                        bed_start = pileup_col.pos
                        curr_flag = 'NONUNIQUE_ALIGNMENTS'

            # END PREVIOUS FLAG
            elif curr_flag:
                _save_bed_line(chrom, bed_start, pileup_col.pos, curr_flag)
                curr_flag = None

            next_pos = pileup_col.pos +1

        # END LAST FLAG
        _save_bed_line(chrom, bed_start, pileup_col.pos-1, curr_flag)

    # combine and sort all bed lines by position
    all_bed_lines = list(itertools.chain(*bed_lines.values()))
    all_bed_lines.sort(key=lambda l: l['start'])

    with open(bed_output, 'w') as fh:
        for line in all_bed_lines:

                assert line['start'] <= line['end']

                print >> fh, '{}\t{}\t{}\t{}'.format(
                    line['chrom'],
                    line['start'],
                    line['end'],
                    line['name'])

if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) >= 4:
        args[3] = int(args[3])
    if len(args) >= 5:
        args[4] = int(args[4])

    get_callable_loci(*args)

