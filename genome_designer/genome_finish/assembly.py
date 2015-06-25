import datetime
import os

from Bio import SeqIO

from genome_finish.millstone_de_novo_fns import add_paired_mates
from genome_finish.millstone_de_novo_fns import get_clipped_reads
from genome_finish.millstone_de_novo_fns import get_unmapped_reads
from genome_finish.millstone_de_novo_fns import get_split_reads
from genome_finish.millstone_de_novo_fns import run_velvet
from main.models import Contig
from main.models import Dataset
from main.model_utils import get_dataset_with_type
from pipeline.read_alignment import align_with_bwa_mem
from utils.bam_utils import concatenate_bams
from utils.bam_utils import make_bam
from utils.bam_utils import make_sam
from utils.bam_utils import rmdup
from utils.import_util import add_dataset_to_entity
from utils.import_util import prepare_ref_genome_related_datasets

# Default args for velvet assembly
VELVET_COVERAGE_CUTOFF = 3
VELVET_KMER_LIST = [21]


def generate_contigs(experiment_sample_to_alignment, contig_label_base):

    timestamp = str(datetime.datetime.now())

    # Grab reference genome fasta path
    reference_genome = (
            experiment_sample_to_alignment.alignment_group.reference_genome)
    ref_fasta_dataset = reference_genome.dataset_set.get_or_create(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA)[0]
    prepare_ref_genome_related_datasets(reference_genome, ref_fasta_dataset)

    # Make data_dir directory to house genome_finishing files
    genome_finishing_dir = os.path.join(
            experiment_sample_to_alignment.get_model_data_dir(),
            'genome_finishing')

    # Make data_dir directory if it does not exist
    if not os.path.exists(genome_finishing_dir):
        os.mkdir(genome_finishing_dir)

    data_dir = os.path.join(genome_finishing_dir, '0')
    data_dir_counter = 0
    while(os.path.exists(data_dir)):
        data_dir_counter += 1
        data_dir = os.path.join(genome_finishing_dir, str(data_dir_counter))
    os.mkdir(data_dir)

    # Retrieve bwa mem .bam alignment if exists otherwise generate it
    if not experiment_sample_to_alignment.dataset_set.filter(
            type=Dataset.TYPE.BWA_ALIGN).exists():
        add_dataset_to_entity(
                experiment_sample_to_alignment,
                'sample_alignment_for_assembly',
                Dataset.TYPE.BWA_ALIGN)
        align_with_bwa_mem(
                experiment_sample_to_alignment.alignment_group,
                experiment_sample_to_alignment,
                project=reference_genome.project)

    # Grab alignment bam file-path
    alignment_sorted = get_dataset_with_type(
                experiment_sample_to_alignment,
                Dataset.TYPE.BWA_ALIGN).get_absolute_location()

    alignment_prefix = os.path.join(data_dir, 'bwa_align')

    # Extract SV indicating reads
    unmapped_reads = alignment_prefix + '.unmapped.bam'
    get_unmapped_reads(alignment_sorted, unmapped_reads)

    split_reads = alignment_prefix + '.split.bam'
    get_split_reads(alignment_sorted, split_reads)

    clipped_reads = alignment_prefix + '.clipped.bam'
    get_clipped_reads(alignment_sorted, clipped_reads)

    # Aggregate SV indicants
    SV_indicants_bam = alignment_prefix + '.SV_indicants.bam'
    concatenate_bams(
            [unmapped_reads, split_reads, clipped_reads],
            SV_indicants_bam)

    # Remove duplicates
    SV_indicants_no_dups_bam = alignment_prefix + '.SV_indicants_no_dups.bam'
    rmdup(SV_indicants_bam, SV_indicants_no_dups_bam)

    # Convert SV indicants bam to sam
    SV_indicants_sam = alignment_prefix + '.SV_indicants.sam'
    make_sam(SV_indicants_no_dups_bam, SV_indicants_sam)

    # Add mate pairs to SV indicants sam
    SV_indicants_with_pairs_sam = (
        alignment_prefix + '.SV_indicants_with_pairs.sam')
    add_paired_mates(
            SV_indicants_sam, alignment_sorted, SV_indicants_with_pairs_sam)

    # Make bam of SV indicants w/mate pairs
    SV_indicants_with_pairs_bam = (
        alignment_prefix + '.SV_indicants_with_pairs.bam')
    make_bam(SV_indicants_with_pairs_sam, SV_indicants_with_pairs_bam)

    # Velvet assembly
    contig_files = []
    opt_dict = {
        'velveth': {
            'shortPaired': ''
        },
        'velvetg': {
            'cov_cutoff': VELVET_COVERAGE_CUTOFF
        }
    }
    kmer_list = VELVET_KMER_LIST
    for kmer_length in kmer_list:
        # Set hash length argument for velveth
        opt_dict['velveth']['hash_length'] = kmer_length

        # Run velvet assembly on SV indicants
        velvet_dir = os.path.join(data_dir, 'velvet_k' + str(kmer_length))
        run_velvet(
                SV_indicants_with_pairs_bam,
                velvet_dir,
                opt_dict)

        # Collect resulting contigs fasta
        contigs_fasta = os.path.join(velvet_dir, 'contigs.fa')
        contig_files.append(contigs_fasta)

        for seq_record in SeqIO.parse(contigs_fasta, 'fasta'):

            contig_label = '_'.join(
                    [contig_label_base, seq_record.description])

            # Create an insertion model for the contig
            contig = Contig.objects.create(
                    label=contig_label,
                    parent_reference_genome=reference_genome,
                    experiment_sample_to_alignment=(
                            experiment_sample_to_alignment))

            contig.metadata['coverage'] = float(
                    seq_record.description.rsplit('_', 1)[1])
            contig.metadata['timestamp'] = timestamp

            contig.ensure_model_data_dir_exists()

            dataset_path = os.path.join(contig.get_model_data_dir(),
                    'fasta.fa')
            assert not os.path.exists(dataset_path)

            with open(dataset_path, 'w') as fh:
                SeqIO.write([seq_record], fh, 'fasta')

            add_dataset_to_entity(
                    contig,
                    'contig_fasta',
                    Dataset.TYPE.REFERENCE_GENOME_FASTA,
                    filesystem_location=dataset_path)

    return contig_files
