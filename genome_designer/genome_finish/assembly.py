import sys
import os

from genome_finish import de_novo_bam
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome

def truncate_ext(path, count=1):
    count = count -1
    if count == 0:
        return path[:path.rindex(".")]
    else:
        return truncate_ext(path[:path.rindex(".")], count)

def generate_contigs(experimentSampleToAlignment):

    # Get fasta read files
    experimentSample = experimentSampleToAlignment.experiment_sample
    fastq1 = experimentSample.dataset_set.get(type=Dataset.TYPE.FASTQ1).get_absolute_location()
    fastq2 = experimentSample.dataset_set.get(type=Dataset.TYPE.FASTQ2).get_absolute_location()
    read_fastqs = [fastq1, fastq2]

    # Get fasta reference genome file
    referenceGenome = experimentSampleToAlignment.alignment_group.reference_genome
    ref_fasta = referenceGenome.dataset_set.get(type = Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    
    # Make data_dir directory to house genome_finishing files
    project_dir = referenceGenome.project.get_model_data_dir()
    data_dir = os.path.join(project_dir, "genome_finishing", experimentSampleToAlignment.uid)

    # Make data_dir directory if it does not exist
    if not os.path.exists(os.path.join(project_dir, "genome_finishing")):
        os.mkdir(os.path.join(project_dir, "genome_finishing"))
        os.mkdir(data_dir)
    elif not os.path.exists(data_dir):
        os.mkdir(data_dir)


    # ---------- Begin Genome Finishing ---------- #

    # Align reads
    alignment_sorted = de_novo_bam.bwa_align(read_fastqs, ref_fasta, data_dir)

    # Trim 'alignment.sorted.bam' off filename to get bam_prefix
    bam_prefix = truncate_ext(alignment_sorted, 3)

    # Get split and unmapped reads
    unmapped_reads = bam_prefix + ".alignment.unmapped.bam"
    # alignment_sorted = bam_prefix + ".alignment.sorted.bam"
    de_novo_bam.millstone_de_novo_fns.get_unmapped_reads(bam_prefix + ".alignment.sorted.bam", unmapped_reads)

    split_reads = bam_prefix + ".alignment.split.bam"
    de_novo_bam.millstone_de_novo_fns.get_split_reads(bam_prefix + ".alignment.sorted.bam", split_reads)

    # Concatenate split and unmapped reads
    unmapped_and_split = bam_prefix + ".alignment.unmapped_and_split.bam"
    de_novo_bam.concatenate_bams(unmapped_reads, split_reads, unmapped_and_split)

    #TODO: sort?

    print "Make split_unmapped sam"
    # Make split_unmapped sam (for add_paired_mates)
    unmapped_and_split_sam = bam_prefix + ".alignment.unmapped_and_split.sam"
    de_novo_bam.make_sam(unmapped_and_split, unmapped_and_split_sam)

    print "Add mate pairs"
    # Add mate pairs %%%%%%% samtools header error if try to view unmapped_and_split.sam, but can view bam
    unmapped_and_split_with_pairs_sam = bam_prefix + ".alignment.unmapped_and_split_with_pairs.sam"
    de_novo_bam.millstone_de_novo_fns.add_paired_mates(unmapped_and_split_sam, alignment_sorted, unmapped_and_split_with_pairs_sam)

    print "Make bam of mate pairs"
    # Make bam of mate pairs %%%%%%% samtools header error
    unmapped_and_split_with_pairs_bam = bam_prefix + ".alignment.unmapped_and_split_with_pairs.bam"
    de_novo_bam.make_bam(unmapped_and_split_with_pairs_sam, unmapped_and_split_with_pairs_bam)

    #TODO: sort?

    print "Remove duplicates"
    # Remove duplicates %%%%%% seems to delete too much
    #de_novo_bam.rmdup(unmapped_and_split_with_pairs_bam) 

    print "Assemble with velvet"
    # Assemble with velvet

    # contig_alignment_dir = os.path.join(data_dir, "contig_alignments")
    # os.mkdir(contig_alignment_dir)

    kmer_list = [21] #range(11,33,2)
    contig_files = []

    for kmer_length in kmer_list:
        print "assembling with kmer length " + str(kmer_length)
        velvet_dir = os.path.join(data_dir, "velvet_k" + str(kmer_length))
        de_novo_bam.run_velvet(
                unmapped_and_split_with_pairs_bam,
                output_dir_name = velvet_dir,
                hash_length = kmer_length,
                no_defaults = True) #Change back to False when fastq's more realistic

        contigs_fasta = os.path.join(velvet_dir,"contigs.fa")
        contig_files.append(contigs_fasta)

    return contig_files

        # # Align contigs
        # contigs = os.path.join(velvet_dir, "contigs.fa")
        # de_novo_bam.bwa_align(
        #         [contigs],
        #         actual_genome,
        #         os.path.join(contig_alignment_dir, "k" + str(kmer_length) + ".alignment.bam")
        # )
