"""Functions for attempting to de novo align unmapped reads.
"""

import os
import shutil
import subprocess

from django.conf import settings

from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.models import Dataset
from main.models import get_or_create_derived_bam_dataset
from pipeline.read_alignment import get_split_reads


VELVETH_BINARY = '%s/velvet/velveth' % settings.TOOLS_DIR
VELVETG_BINARY = '%s/velvet/velvetg' % settings.TOOLS_DIR
VELVET_HASH_LENGTH = 21


def get_unmapped_reads(sample_alignment, force_rerun=False):
    """Returns Dataset for unmapped reads from the sample alignment.

    Computes it if it doesn't exist already.

    Args:
        sample_alignment: Complete sample_alignment.
        force_rerun: If True, Dataset is recomputed.

    Returns:
        Dataset that points to the bam file containing unmapped reads from
        sample_alignment.
    """
    def derivation_fn(sample_alignment, unmapped_reads_dataset):
        # Get the original bam file.
        bam_dataset = get_dataset_with_type(sample_alignment,
                Dataset.TYPE.BWA_ALIGN)
        bam_filename = bam_dataset.get_absolute_location()

        # Allocate a filename for the unmapped reads.
        unmapped_reads_bam_file = (os.path.splitext(bam_filename)[0] +
                '.unmapped.bam')
        unmapped_reads_dataset.filesystem_location = clean_filesystem_location(
                unmapped_reads_bam_file)
        unmapped_reads_dataset.save(update_fields=['filesystem_location'])

        cmd = '{samtools} view -h -b -f 0x4 {bam_filename}'.format(
                samtools=settings.SAMTOOLS_BINARY,
                bam_filename=bam_filename)
        with open(unmapped_reads_bam_file, 'w') as output_fh:
           subprocess.check_call(cmd, stdout=output_fh, shell=True)

    return get_or_create_derived_bam_dataset(sample_alignment,
            Dataset.TYPE.BWA_UNMAPPED, derivation_fn, force_rerun=force_rerun)


def get_bam_for_de_novo_alignment(sample_alignment,
        force_include_reads_in_intervals=[],
        force_rerun=False):
    """Prepares a bam file for de novo alignment.

    Args:
        sample_alignment: Complete sample_alignment.
        force_include_reads_in_intervals: List of intervals where reads near
            those intervals are included.
        force_rerun: If True, Dataset is recomputed.

    Returns:
        Dataset that points to the bam file ready for de novo alignment.
    """
    def derivation_fn(sample_alignment, new_dataset):
        """Creates a bam file with reads relevant for de novo assembly.
        """
        # Get the original bam file.
        bam_dataset = get_dataset_with_type(sample_alignment,
                Dataset.TYPE.BWA_ALIGN)
        orig_bam_filename = bam_dataset.get_absolute_location()

        # Allocate a filename for the new dataset.
        de_novo_bam_filelocation = (os.path.splitext(orig_bam_filename)[0] +
                '.de_novo.bam')
        new_dataset.filesystem_location = clean_filesystem_location(
                de_novo_bam_filelocation)
        new_dataset.save(update_fields=['filesystem_location'])

        ### Strategy
        # 0. Create intermediate sam
        # 1. Grab unmapped reads.
        #     a. If the corresponding pair for any read was mapped (and thus not
        #        in the unmapped dataset), grab that from the original bam file.
        # 2. Grab split reads.
        # 2b. Add reads that are in intervals of interest.
        # 3. Sort by name so that paired reads are next to each other
        # 4. Filter out duplicate reads. Requires sort in step 3.
        # 5. Convert to bam.
        # 6. Delete intermediate files.

        # 0. Create intermediate files.
        intermediate_sam = (
                os.path.splitext(de_novo_bam_filelocation)[0] + '.int.sam')
        paired_intermediate_sam = (
                os.path.splitext(intermediate_sam)[0] + '.paired.sam')
        sorted_intermediate_sam = (
                os.path.splitext(paired_intermediate_sam)[0] + '.sorted.sam')
        deduped_sorted_intermediate_sam = (
                os.path.splitext(sorted_intermediate_sam)[0] + '.deduped.sam')
        intermediate_files = [
                intermediate_sam,
                paired_intermediate_sam,
                sorted_intermediate_sam,
                deduped_sorted_intermediate_sam
        ]

        # 1. Get unmapped reads.
        unmapped_reads_dataset = get_unmapped_reads(sample_alignment)
        unmapped_reads_bam_file = unmapped_reads_dataset.get_absolute_location()
        cmd = '{samtools} view -h {bam_filename}'.format(
                samtools=settings.SAMTOOLS_BINARY,
                bam_filename=unmapped_reads_bam_file)
        with open(intermediate_sam, 'w') as output_fh:
           subprocess.check_call(cmd, stdout=output_fh, shell=True)

        # 2. Grab split reads.
        split_reads_dataset = get_split_reads(sample_alignment)
        split_rads_bam_file = split_reads_dataset.get_absolute_location()
        cmd = '{samtools} view {bam_filename}'.format(
                samtools=settings.SAMTOOLS_BINARY,
                bam_filename=split_rads_bam_file)
        with open(intermediate_sam, 'a') as output_fh:
           subprocess.check_call(cmd, stdout=output_fh, shell=True)

        # 2b. Add reads that are in any of intervals that we want to include.
        if force_include_reads_in_intervals:
            with open(intermediate_sam, 'a') as output_fh:
                get_reads_in_interval_list(orig_bam_filename,
                        force_include_reads_in_intervals, output_fh)

        # 2c. For each of the reads that we've included, grab their
        # corresponding pairs.
        add_paired_mates(intermediate_sam, orig_bam_filename,
                paired_intermediate_sam)

        # 3. Sort by name so that paired reads are next to each other
        cmd = (
                '(grep ^"@" {sam_file}; grep -v ^"@" {sam_file} | '
                'sort -k1,1 -k2,2n) > {sorted_sam_file}'
        ).format(
                sam_file=paired_intermediate_sam,
                sorted_sam_file=sorted_intermediate_sam)
        subprocess.call(cmd, shell=True)

        # 4. Filter out duplicate reads. Requires sort in step 3.
        cmd = 'uniq {sorted_sam_file} > {deduped_sam_file}'.format(
                sorted_sam_file=sorted_intermediate_sam,
                deduped_sam_file=deduped_sorted_intermediate_sam)
        subprocess.call(cmd, shell=True)

        # 5. Convert to bam.
        cmd = '{samtools} view -Sb {sam_file} > {final_bam_file}'.format(
                samtools=settings.SAMTOOLS_BINARY,
                sam_file=intermediate_sam,
                final_bam_file=de_novo_bam_filelocation)
        subprocess.call(cmd, shell=True)

        # 6. Delete intermediate files.
        for f in intermediate_files:
            os.remove(f)

    return get_or_create_derived_bam_dataset(sample_alignment,
            Dataset.TYPE.BWA_FOR_DE_NOVO_ASSEMBLY, derivation_fn,
            force_rerun=force_rerun)


def get_reads_in_interval_list(bam_filename, interval_list, output_fh):
    """Update output_fh with reads that map to positions in the
    interval_list. Assumes output_fh is a sam file.
    """
    cmd = '{samtools} view {bam_filename}'.format(
            samtools=settings.SAMTOOLS_BINARY,
            bam_filename=bam_filename)
    read_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    for line in read_proc.stdout:
        parts = line.split()
        pos = int(parts[3])
        for interval in interval_list:
            if interval[0] <= pos <= interval[1]:
                output_fh.write(line)
                break


def add_paired_mates(input_sam_path, source_bam_filename, output_sam_path):
    """Creates a file at output_sam_path that contains all the reads in
    input_sam_path, as well as all of their paired mates.

    The resulting sam is not sorted and may contain duplicates. Clients
    should filter elsewhere.

    TODO: This is currently nasty inefficient. Is there a better way to do
    this, e.g. leverage indexing somehow?
    """
    # First, copy the input file to output in order to preserve the header.
    shutil.copyfile(input_sam_path, output_sam_path)

    # First count id occurrences to avoid looking up those that happen more than
    # once.
    cmd = [
            settings.SAMTOOLS_BINARY,
            'view',
            '-S',
            input_sam_path
    ]
    initial_read_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    # Input might have duplicate lines (same element of each pair repeated), so
    # we leverage the fact that flags are unique for each of the elements in a
    # pair to do counting.
    read_id_to_flags_map = {}
    for line in initial_read_proc.stdout:
        parts = line.split()
        read_id = parts[0]
        flags_value = parts[1]
        if not read_id in read_id_to_flags_map:
            read_id_to_flags_map[read_id] = []
        if not read_id in read_id_to_flags_map:
            read_id_to_flags_map[read_id] = flags_value
        assert len(read_id_to_flags_map[read_id]) <= 2

    # Next, read through the output, copying each line, and additionally
    # searching through the bam file for the corresponding paired read.
    # Next, read through the sam input and append lines to the output that
    # are the corresponding mate pairs.
    with open(output_sam_path, 'a') as output_fh:
        cmd = [
                settings.SAMTOOLS_BINARY,
                'view',
                '-S',
                input_sam_path
        ]
        read_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in read_proc.stdout:
            parts = line.split()
            read_id = parts[0]

            # No need to look up if we already have both.
            # if len(read_id_to_flags_map[read_id]) == 2:
            #     continue

            # Find the corresponding reads in the bam file.
            samtools_cmd = [
                    settings.SAMTOOLS_BINARY,
                    'view',
                    source_bam_filename
            ]
            samtools_proc = subprocess.Popen(samtools_cmd,
                    stdout=subprocess.PIPE)
            grep_cmd = [
                    'grep',
                    '-m', '2',
                    read_id
            ]
            grep_proc = subprocess.Popen(grep_cmd,
                    stdout=subprocess.PIPE, stdin=samtools_proc.stdout)
            samtools_proc.stdout.close()
            for bam_line in grep_proc.stdout:
                if bam_line != line:
                    output_fh.write(bam_line)
            grep_proc.terminate()
            samtools_proc.terminate()


def run_velvet(sample_alignment, output_dir_name='velvet',
        force_rerun=False):
    velvet_output_dir = os.path.join(
            sample_alignment.get_model_data_dir(), output_dir_name)

    bam_dataset = get_bam_for_de_novo_alignment(sample_alignment,
            force_rerun=force_rerun)
    bam_file = bam_dataset.get_absolute_location()

    # First run velveth to build hash of reads.
    cmd = '{velveth} {output_dir} {hash_length} -bam {bam_file}'.format(
            velveth=VELVETH_BINARY,
            output_dir=velvet_output_dir,
            hash_length=VELVET_HASH_LENGTH,
            bam_file=bam_file
    )
    subprocess.call(cmd, shell=True)

    # Then run velvetg to build graph and find contigs.
    cmd = '{velvetg} {output_dir} -cov_cutoff 20 -ins_length 200 -ins_length_sd 90'.format(
            velvetg=VELVETG_BINARY,
            output_dir=velvet_output_dir
    )
    subprocess.call(cmd, shell=True)
