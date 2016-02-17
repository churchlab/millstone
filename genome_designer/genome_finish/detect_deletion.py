import numpy as np
import pysam

bam_path = '/home/wahern/projects/millstone/genome_designer/conf/../temp_data/projects/cac10232/alignment_groups/9a10bebe/sample_alignments/59ec2008/bwa_align.sorted.withmd.bam'


def per_base_cov_stats(bam_path):

    bamfile = pysam.AlignmentFile(bam_path, 'rb')

    chrom_list = bamfile.references
    chrom_lens = bamfile.lengths
    c_starts = [0] * len(chrom_list)
    c_ends = chrom_lens

    chrom_to_cov_list = {}

    for chrom, c_start, c_end in zip(chrom_list, c_starts, c_ends):
        altalign_arr = np.zeros(c_end)
        depth_arr = np.zeros(c_end)

        for pileup_col in bamfile.pileup(chrom,
                start=c_start, end=c_end, truncate=True):

            # number of segments aligned to this position
            depth = pileup_col.nsegments
            # number of segments who also map somewhere else
            altaligns = 0

            # fill in badmapq & altaligns
            for p in pileup_col.pileups:
                a = p.alignment
                tag_as = a.get_tag('AS')
                tag_xs = a.get_tag('XS')
                altaligns += tag_as <= tag_xs

            assert altalign_arr[pileup_col.reference_pos] == 0
            assert depth_arr[pileup_col.reference_pos] == 0
            altalign_arr[pileup_col.reference_pos] = altaligns
            depth_arr[pileup_col.reference_pos] = depth

        chrom_to_cov_list[chrom] = {
                'altaligns': altalign_arr,
                'depths': depth_arr
        }

    return chrom_to_cov_list
