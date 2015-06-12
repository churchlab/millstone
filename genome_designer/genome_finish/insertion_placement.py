import os
import tempfile

from Bio import SeqIO
from django.conf import settings

from genome_finish.millstone_de_novo_fns import get_insertion_location
from genome_finish.millstone_de_novo_fns import get_local_contig_placement
from genome_finish.millstone_de_novo_fns import make_sliced_fasta
from main.exceptions import ValidationException
from main.model_utils import get_dataset_with_type
from main.models import Dataset
from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from pipeline.read_alignment import align_with_bwa_mem
from utils import convert_seqrecord_to_fastq
from utils import generate_safe_filename_prefix_from_label
from utils.import_util import add_dataset_to_entity
from utils.reference_genome_maker_util import generate_new_reference_genome


def place_contig(reference_genome, contig_seqrecord,
        new_reference_genome_label):
    # Find the insertion position in reference and the end positions of the
    # insertion cassette in the contig
    placement_position_params = find_contig_insertion_site(
            reference_genome,
            contig_seqrecord)

    new_reference_genome_params = {
        'label': new_reference_genome_label
    }

    # Generate a new version of the reference genome with the
    # cassette incorporated
    return place_cassette(
            reference_genome, contig_seqrecord,
            placement_position_params,
            new_reference_genome_params)


def align_contig_to_reference(reference_genome, contig_seqrecord):
    alignment_group = AlignmentGroup.objects.create(
            reference_genome=reference_genome,
            label='contig_alignment')

    contig_sample = ExperimentSample.objects.create(
            project=reference_genome.project,
            label='contig_sample_2')

    # Convert inserted sequence to fastq for alignment
    fastq_path = os.path.join(
            contig_sample.get_model_data_dir(), 'insertion.fq')
    convert_seqrecord_to_fastq(contig_seqrecord, fastq_path)

    # Add fastq dataset to Experiment Sample
    add_dataset_to_entity(
            contig_sample, 'contig_fastq', Dataset.TYPE.FASTQ1,
            filesystem_location=fastq_path)

    sample_to_alignment = ExperimentSampleToAlignment.objects.create(
            alignment_group=alignment_group,
            experiment_sample=contig_sample)

    add_dataset_to_entity(
            sample_to_alignment, 'contig_to_ref_bam', Dataset.TYPE.BWA_ALIGN)

    align_with_bwa_mem(
            alignment_group,
            sample_to_alignment,
            project=reference_genome.project)

    contig_to_ref_bam = get_dataset_with_type(
            sample_to_alignment, Dataset.TYPE.BWA_ALIGN
                    ).get_absolute_location()

    return contig_to_ref_bam


def find_contig_insertion_site(reference_genome, contig_seqrecord):

    # Align contig to reference
    contig_to_ref_bam = align_contig_to_reference(
            reference_genome, contig_seqrecord)

    # Find region of insertion in the reference genome
    insertion_location_data = get_insertion_location(contig_to_ref_bam)

    ref_chromosome_seqrecord_id = insertion_location_data['chromosome_seqrecord_id']
    ins_left_end = insertion_location_data['left_end']
    ins_right_end = insertion_location_data['right_end']

    # Make temporary FASTA for only the region of insertion in the reference
    reference_genome_fasta = get_dataset_with_type(
            reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    ref_slice_fasta_temp = _make_temp_file(
            '_'.join([reference_genome.label, 'slice']), '.fa')

    make_sliced_fasta(reference_genome_fasta,
            ref_chromosome_seqrecord_id, ins_left_end, ins_right_end,
            ref_slice_fasta_temp)

    # Make temporary FASTA for contig_seqrecord
    contig_fasta_temp = _make_temp_file(contig_seqrecord.id, '.fa')
    SeqIO.write(contig_seqrecord, contig_fasta_temp, 'fasta')

    # Use Seqan C++ script to find exact insertion position in reference
    # and the exact ends of the insertion cassette in the contig
    local_contig_placement = get_local_contig_placement(ref_slice_fasta_temp,
            contig_fasta_temp)

    # The insertion position in the reference genome is the position of the
    # cut out region of insertion and the exact offset within it found by
    # the split alignment Seqan script
    ref_genome_insertion_pos = (
            ins_left_end + local_contig_placement['ref_split_pos'] + 1)

    return {
        'ref_insertion_pos': ref_genome_insertion_pos,
        'ref_chromosome_seqrecord_id': ref_chromosome_seqrecord_id,
        'contig_cassette_start_pos': local_contig_placement['contig_start_pos'],
        'contig_cassette_end_pos': local_contig_placement['contig_end_pos']
    }


def _make_temp_file(label, extension):
    if not os.path.exists(settings.TEMP_FILE_ROOT):
        os.mkdir(settings.TEMP_FILE_ROOT)
    _, temp_file_path = tempfile.mkstemp(
            suffix=('_' + generate_safe_filename_prefix_from_label(label) +
                    extension),
            dir=settings.TEMP_FILE_ROOT)
    return temp_file_path


def place_cassette(reference_genome, contig_seqrecord,
        placement_position_params, new_reference_genome_params):

    # Validate param dictionaries
    try:
        for key in ['ref_insertion_pos', 'ref_chromosome_seqrecord_id',
                'contig_cassette_start_pos', 'contig_cassette_end_pos']:
            assert key in placement_position_params

        assert 'label' in new_reference_genome_params
    except Exception as e:
        raise ValidationException(e)

    # Get chromosome of reference genome to be recieving cassette
    insertion_chromosome = Chromosome.objects.get(
        reference_genome=reference_genome,
        seqrecord_id=placement_position_params['ref_chromosome_seqrecord_id'])

    # Create variant to house insertion
    insertion_variant = Variant.objects.create(
            reference_genome=reference_genome,
            chromosome=insertion_chromosome,
            type=Variant.TYPE.INSERTION,
            position=placement_position_params['ref_insertion_pos'],
            ref_value='')

    # Extract cassette sequence from contig
    cassette_sequence = str(contig_seqrecord.seq[
            placement_position_params['contig_cassette_start_pos']:
            placement_position_params['contig_cassette_end_pos']])

    # Create the variant alternate for the cassette sequence
    insertion_variant.variantalternate_set.add(
            VariantAlternate.objects.create(
                    variant=insertion_variant,
                    alt_value=cassette_sequence))

    # House insertion variant in variant set to be applied to the ref
    insertion_variant_set = VariantSet.objects.create(
            reference_genome=reference_genome,
            label='insertion_variant_set')

    VariantToVariantSet.objects.create(
            variant=insertion_variant,
            variant_set=insertion_variant_set)

    return generate_new_reference_genome(
            insertion_variant_set,
            new_reference_genome_params)
