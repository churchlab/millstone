"""
Utility functions for working with bam files.
"""

import os
import shutil
import subprocess

from django.conf import settings
import pysam
import numpy as np

from utils import convert_fasta_to_fastq

BWA_BINARY = os.path.join(settings.TOOLS_DIR, 'bwa/bwa')


def clipping_stats(bam_path, sample_size=1000):

    BAM_CSOFT_CLIP = 4
    BAM_CHARD_CLIP = 5
    clip_codes = [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]

    samfile = pysam.AlignmentFile(bam_path)
    sample_size = min(sample_size, samfile.mapped)

    terminal_clipping = []
    i = 0
    for read in samfile:
        if not read.is_unmapped:
            first_cig = read.cigartuples[0]
            last_cig = read.cigartuples[-1]
            terminal_clipping.append(max([
                    first_cig[1] if first_cig[0] in clip_codes else 0,
                    last_cig[1] if last_cig[0] in clip_codes else 0]))

            i += 1
            if i == sample_size:
                break

    return {'mean': np.mean(terminal_clipping),
            'std': np.std(terminal_clipping)}


def index_bam(bam):
    cmd = "{samtools} index {bam}".format(
            samtools=settings.SAMTOOLS_BINARY,
            bam=bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def sort_bam_by_name(input_bam, output_bam=None):
    if output_bam is None:
        output_bam = input_bam

    cmd = "{samtools} sort -n {input_bam} {output_bam_prefix}".format(
            samtools=settings.SAMTOOLS_BINARY,
            input_bam=input_bam,
            output_bam_prefix=os.path.splitext(output_bam)[0])

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def sort_bam_by_coordinate(input_bam, output_bam=None):
    if output_bam is None:
        output_bam = input_bam

    cmd = "{samtools} sort {input_bam} {output_bam_prefix}".format(
            samtools=settings.SAMTOOLS_BINARY,
            input_bam=input_bam,
            output_bam_prefix=os.path.splitext(output_bam)[0])

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def make_sam(bam, sam_filename=None):
    if sam_filename is None:
        sam_filename = os.path.splitext(bam)[0] + ".sam"

    cmd = "{samtools} view -h {bam} > {sam}".format(
            samtools=settings.SAMTOOLS_BINARY,
            bam=bam,
            sam=sam_filename)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def make_bam(sam, bam_filename=None):
    if bam_filename is None:
        bam_filename = os.path.splitext(sam)[0] + ".bam"

    cmd = "{samtools} view -b -S {sam} > {bam}".format(
            samtools=settings.SAMTOOLS_BINARY,
            sam=sam,
            bam=bam_filename)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def concatenate_bams(bam_list, output):

    cmd = "{samtools} cat -o {output} {bam_files}".format(
            samtools=settings.SAMTOOLS_BINARY,
            bam_files=" ".join(bam_list),
            output=output)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def rmdup(input_bam_file, output_bam_file):
    """ Remove duplicate lines from a SAM file, implemented with sort/uniq
        as samtools 0.1.20 has known bugs with its rmdup function
    """
    # Store input sam header
    output_sam = os.path.splitext(output_bam_file)[0] + ".sam"
    subprocess.check_call(
            ' '.join(
                    [settings.SAMTOOLS_BINARY, 'view', '-H', '-o',
                     output_sam,  input_bam_file]),
            shell=True, executable=settings.BASH_PATH)

    # Convert input bam to sam, sort, remove duplicate adjacent lines,
    # and append to header
    cmd = ' | '.join([
            settings.SAMTOOLS_BINARY + ' view ' + input_bam_file,
            'sort',
            'uniq'
            ]) + ' >> ' + output_sam
    subprocess.check_call(cmd, shell=True, executable=settings.BASH_PATH)

    make_bam(output_sam, output_bam_file)


def filter_bam_file_by_row(input_bam_path, filter_fn, output_bam_path):
    """Filters rows out of a bam file that don't pass a given filter function.

    This function keeps all header lines.

    Args:
        input_bam_path: Absolute path to input bam file.
        filter_fn: Function applied to each row of the input bam and returns a
            Boolean. If True, keeps the row.
        output_bam_path: Absolute path to the output bam file.
    """
    output_root = os.path.splitext(output_bam_path)[0]
    initial_sam_intermediate = output_root + '.sam'
    filtered_sam_intermediate = output_root + '.filtered.sam'
    final_bam = output_root + '.filtered.bam'

    # Convert to SAM (preserve header with -h option).
    with open(initial_sam_intermediate, 'w') as output_fh:
        subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-h', input_bam_path],
                stdout=output_fh)

    # Filter.
    with open(filtered_sam_intermediate, 'w') as output_fh:
        with open(initial_sam_intermediate) as input_fh:
            for line in input_fh:
                # Always write header lines.
                if line[0] == '@':
                    output_fh.write(line)
                    continue

                if filter_fn(line):
                    output_fh.write(line)
                    continue

    # Write final bam.
    with open(final_bam, 'w') as fh:
        subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-bS',
                 filtered_sam_intermediate],
                stdout=fh)

    # Move temp file to the original file location.
    shutil.move(final_bam, output_bam_path)

    # Delete intermediate files.
    os.remove(initial_sam_intermediate)
    os.remove(filtered_sam_intermediate)


def minimal_bwa_align(reads, ref_fasta, data_dir):
    # 0. Interpret reads file type
    if all([r.endswith(".fa") for r in reads]):
        reads_fq = [os.path.join(data_dir, os.path.splitext(r)[0] + ".fq")
                    for r in reads]
        for i, r in enumerate(reads):
            convert_fasta_to_fastq(r, reads_fq[i])
    elif all([r.endswith(".fq") for r in reads]):
        reads_fq = reads
    else:
        raise(Exception("All reads must have file extension .fq or .fa"))

    filename_prefix = os.path.join(data_dir, "bwa_align.alignment")

    # 1. bwa index ref.fa #TODO: Check if already indexed
    cmd = "{bwa} index {ref_path}".format(
            bwa=BWA_BINARY,
            ref_path=ref_fasta)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 2. bwa mem ref.fa contigs.fq > alignment.sam
    alignment_sam = filename_prefix + ".sam"

    cmd = "{bwa} mem {ref_fasta} {contigs_fastq} > {alignment_sam}".format(
            bwa=BWA_BINARY,
            contigs_fastq=" ".join(reads_fq),
            ref_fasta=ref_fasta,
            alignment_sam=alignment_sam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 3. samtools view -b alignment.sam > alignment.bam
    alignment_bam = filename_prefix + ".bam"

    cmd = "{samtools} view -b -S {alignment_sam} > {alignment_bam}".format(
            samtools=settings.SAMTOOLS_BINARY,
            alignment_sam=alignment_sam,
            alignment_bam=alignment_bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 4. samtools sort alignment.bam
    sort_bam_by_coordinate(alignment_bam)

    # 5. Index it
    index_bam(alignment_bam)

    return alignment_bam
