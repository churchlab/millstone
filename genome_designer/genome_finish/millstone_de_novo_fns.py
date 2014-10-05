import subprocess
import os
import shutil
import re

from settings import TOOLS_DIR
from settings import SAMTOOLS_BINARY
from settings import BASH_PATH

def get_unmapped_reads(bam_filename, output_filename):

    unmapped_reads_bam_file = output_filename

    cmd = '{samtools} view -h -b -f 0x4 {bam_filename}'.format(
            samtools=SAMTOOLS_BINARY,
            bam_filename=bam_filename)
    with open(output_filename, 'w') as output_fh:
       subprocess.call(cmd, stdout=output_fh, shell=True, executable=BASH_PATH)


def get_split_reads(bam_filename, output_filename):
    """Isolate split reads from a sample alignment.

    This uses a python script supplied with Lumppy, that is run as a
    separate process.

    NOTE THAT THIS SCRIPT ONLY WORKS WITH BWA MEM.
    """

    # Use lumpy bwa-mem split read script to pull out split reads.
    filter_split_reads = ' | '.join([
            '{samtools} view -h {bam_filename}',
            'python {lumpy_bwa_mem_sr_script} -i stdin',
            '{samtools} view -Sb -']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    lumpy_bwa_mem_sr_script= os.path.join(
                            TOOLS_DIR, 'lumpy','extractSplitReads_BwaMem'))

    try:
		fh = open(output_filename, 'w')
		subprocess.check_call(filter_split_reads, stdout=fh, shell=True, executable=BASH_PATH)
		fh.close()

        # sort the split reads, overwrite the old file
		subprocess.check_call([SAMTOOLS_BINARY, 'sort', output_filename,
                os.path.splitext(output_filename)[0]])

    except subprocess.CalledProcessError:
		raise Exception("Exception caught in split reads generator, perhaps due to no split reads")


def _parse_sam_line(line):
    parts = line.split()
    return {
        'read_id': parts[0],
        'flags': parts[1]
    }


def add_paired_mates(input_sam_path, source_bam_filename, output_sam_path):
    """Creates a file at output_sam_path that contains all the reads in
    input_sam_path, as well as all of their paired mates.

    The resulting sam is not sorted and may contain duplicates. Clients
    should filter elsewhere.

    TODO: This is currently nasty inefficient. Is there a better way to do
    this, e.g. leverage indexing somehow?

    TODO: This is potentially memory-overwhelming. Need to think about
    the de novo assembly feature more carefully before making it user-facing.
    """
    ### Strategy:
    # 1. Copy input to output, to preserve existing reads.
    # 2. Load all ids into a dictionary, storing information about whether
    #    both pairs are already in the output.
    #    NOTE: This step could potentially overwhelm memory for large datasets.
    # 3. Loop through the bam file, checking against the id dictionary to
    #    see whether or not we want to keep that data.

    # 1. Copy input to output, to preserve existing reads.
    shutil.copyfile(input_sam_path, output_sam_path)

    # 2. Create dictionary of read ids.
    read_id_to_flags_map = {}
    with open(input_sam_path) as fh:
        for line in fh:
            #print "line to parse:   ", line
            if re.match('@', line):
                continue
            #print "line did not match re '@'"
            parsed_line = _parse_sam_line(line)
            read_id = parsed_line['read_id']
            flags_value = parsed_line['flags']
            if not read_id in read_id_to_flags_map:
                read_id_to_flags_map[read_id] = []
            if not flags_value in read_id_to_flags_map[read_id]:
                read_id_to_flags_map[read_id].append(flags_value)
            assert 1 <= len(read_id_to_flags_map[read_id]) <= 2

    # 3. Loop through bam file, appending lines whose pairs aren't already
    #    present.

    with open(output_sam_path, 'a') as output_fh:
        # Use samtools to read the bam.
        read_bam_cmd = [
                SAMTOOLS_BINARY, # changed from settings.SAMTOOLS_BINARY
                'view',
                source_bam_filename
        ]
        samtools_proc = subprocess.Popen(read_bam_cmd, stdout=subprocess.PIPE)

        for bam_line in samtools_proc.stdout:
            # Use flags as unique identifier for each read (i.e. which of the
            # pairs for a read id that we are looking at.
            parsed_line = _parse_sam_line(bam_line)
            read_id = parsed_line['read_id']
            flags = parsed_line['flags']
            if not read_id in read_id_to_flags_map:
                continue
            observed_flags = read_id_to_flags_map[read_id]
            if not flags in observed_flags:
                output_fh.write(bam_line)
                # Update flags in case we encounter same one again.
                observed_flags.append(flags)
