"""
Utility methods for running simLibrary and simNGS-type stuff.
"""

import os
import re
import subprocess


def run_simLibrary(source_fasta, output_library_fasta):
    """Runs sim library on the source_fasta and writig to the output file."""
    # We store the location of the simNGS binaries in the django settings.
    import settings

    source_fasta_fh = open(source_fasta)
    output_fh = open(output_library_fasta, 'w')

    simLibrary_binary = os.path.join(settings.SIM_NGS_BIN, 'simLibrary')
    subprocess.call([
        simLibrary_binary,
        '--paired',
        '--readlen', '500',
        '--coverage', '50',
    ], stdin=source_fasta_fh, stdout=output_fh)

    source_fasta_fh.close()
    output_fh.close()


def run_paired_simNGS(source_library_fasta, output_fq):
    """Runs simNGS to produce fake reads.

    Returns the pair of output fq files, one for each read direction.
    """
    # We store the location of the simNGS binaries in the django settings.
    import settings

    source_fasta_fh = open(source_library_fasta)
    output_fh = open(output_fq, 'w')

    simNGS_binary = os.path.join(settings.SIM_NGS_BIN, 'simNGS')
    subprocess.call([
        simNGS_binary,
        '-p', 'paired',
        settings.SIM_NGS_NOISE_SOURCE
    ], stdin=source_fasta_fh, stdout=output_fh)

    source_fasta_fh.close()
    output_fh.close()

    # Split the output_fq according to read polarity and write to two files.
    output_fh = open(output_fq)
    output_fq_1 = os.path.splitext(output_fq)[0] + '.1.fq'
    output_fq_2 = os.path.splitext(output_fq)[0] + '.2.fq'
    output_fq_1_fh = open(output_fq_1, 'w')
    output_fq_2_fh = open(output_fq_2, 'w')

    # The .fq file format has groups of four lines, and the presence of either
    # /1 or /2 in the first line.
    all_fq_lines = output_fh.readlines()
    output_fh.close()
    cursor = 0
    while cursor < len(all_fq_lines):
        CHUNK_SIZE = 4
        # Read 4 lines.
        current_chunk = all_fq_lines[cursor:cursor + CHUNK_SIZE]

        # Determine polarity from first line.
        if re.match(r'@Frag.*/1', current_chunk[0]):
            write_to = output_fq_1_fh
        elif re.match(r'@Frag.*/2', current_chunk[0]):
            write_to = output_fq_2_fh
        else:
            raise ValueError("Unexpected chunk: " + ''.join(current_chunk))
        write_to.write(''.join(current_chunk))
        cursor += CHUNK_SIZE

    output_fq_1_fh.close()
    output_fq_2_fh.close()

    return [output_fq_1, output_fq_2]
