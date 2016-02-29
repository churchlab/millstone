"""Script to combine ReferenceGenomes into a single one.
"""

import os

from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna

from main.models import Chromosome
from main.models import Dataset
from main.models import ReferenceGenome
from utils.import_util import add_dataset_to_entity
from utils import generate_safe_filename_prefix_from_label
from utils import remove_whitespace

DATASET_TO_SEQIO_FORMAT = {
    Dataset.TYPE.REFERENCE_GENOME_GENBANK: 'genbank',
    Dataset.TYPE.REFERENCE_GENOME_FASTA: 'fasta'
}


def combine_list_allformats(
        reference_genome_list, new_ref_genome_label, project):
    """Combine ReferenceGenomes into a new single ReferenceGenome
    composed of the component parts.

    Args:
        reference_genome_list: List of ReferenceGenome objects.
        new_ref_genome_label: Label for the new ReferenceGenome.
        project: Project to which the new ReferenceGenome will be added.

    Returns:
        Object with keys:
            * is_success
            * new_reference_genome (when is_success = True)
            * error_msg (when is_success = False)
    """
    rg_dataset_list = []
    for ref_genome in reference_genome_list:
        rg_dataset_tup = None
        for dataset_type in [
                Dataset.TYPE.REFERENCE_GENOME_GENBANK,
                Dataset.TYPE.REFERENCE_GENOME_FASTA]:
            filter_result = ref_genome.dataset_set.filter(type=dataset_type)
            if len(filter_result):
                rg_dataset_tup = (ref_genome, filter_result[0])
                break
        if (not rg_dataset_tup or
                not os.path.exists(rg_dataset_tup[1].get_absolute_location())):
            return {
                'is_success': False,
                'error_msg': 'All reference genomes must have an associated \
                        FASTA or Genbank dataset'
            }
        else:
            rg_dataset_list.append(rg_dataset_tup)
    assert len(rg_dataset_list) == len(reference_genome_list)

    # Read the datasets into Biopython SeqRecord objects.
    rg_seqrecord_list = []
    seqrecord_ids = []
    seqrecord_descriptions = []
    for rg, dataset in rg_dataset_list:
        with open(dataset.get_absolute_location()) as input_fh:
            for record in SeqIO.parse(
                    input_fh, DATASET_TO_SEQIO_FORMAT[dataset.type]):
                rg_seqrecord_list.append((rg, record))
                seqrecord_ids.append('_'.join([
                        remove_whitespace(rg.label)[:7],
                        remove_whitespace(record.id)[:8]]))
                seqrecord_descriptions.append(record.description)

    # Create a new ReferenceGenome.
    new_ref_genome = ReferenceGenome.objects.create(
            project=project,
            label=new_ref_genome_label)

    # If ReferenceGenome label and Chromosome id are the same, there will be
    # duplicate seqrecord_ids: resolve by including numeric prefix in id
    seq_record_list = []
    MAX_LOCUS_NAME_LEN = 16
    unique_id_len = len(str(len(seqrecord_ids)))
    label_len = (MAX_LOCUS_NAME_LEN - 2 - unique_id_len) / 2
    for i, seqrecord_id in enumerate(seqrecord_ids):
        rg, seqrecord = rg_seqrecord_list[i]

        if seqrecord_ids.count(seqrecord_id) == 1:
            unique_seqrecord_id = seqrecord_id
        else:
            unique_seqrecord_id = '_'.join([
                str(i),
                remove_whitespace(rg.label)[:label_len],
                remove_whitespace(seqrecord.id)[:label_len]])

        seqrecord.seq.alphabet = ambiguous_dna
        seqrecord.name = unique_seqrecord_id
        seqrecord.id = unique_seqrecord_id

        if seqrecord_descriptions.count(seqrecord.description) > 1:
            seqrecord.description = ' '.join([
                    seqrecord.description,
                    'from Reference Genome:', rg.label])

        seq_record_list.append(seqrecord)
        Chromosome.objects.create(
                reference_genome=new_ref_genome,
                label=seqrecord.id,
                seqrecord_id=seqrecord.id,
                num_bases=len(seqrecord))

    # Generate a filename from the label with non-alphanumeric characters
    # replaced by underscores.
    filename_prefix = generate_safe_filename_prefix_from_label(
            new_ref_genome_label)
    does_list_include_genbank = (
            Dataset.TYPE.REFERENCE_GENOME_GENBANK in
            [rg_dataset_tup[1].type for rg_dataset_tup in rg_dataset_list])

    if does_list_include_genbank:
        filename = filename_prefix + '.gb'
    else:
        filename = filename_prefix + '.fa'
    new_file_dest = os.path.join(new_ref_genome.get_model_data_dir(), filename)

    # Write the result.
    if does_list_include_genbank:
        ref_genome_dataset_type = Dataset.TYPE.REFERENCE_GENOME_GENBANK
    else:
        ref_genome_dataset_type = Dataset.TYPE.REFERENCE_GENOME_FASTA
    output_file_format = DATASET_TO_SEQIO_FORMAT[ref_genome_dataset_type]

    with open(new_file_dest, 'w') as output_fh:
        SeqIO.write(seq_record_list, output_fh, output_file_format)

    # Create a dataset which will point to the file. This step must happen
    # after writing the file because a signal will be triggered which requires
    # the Genbank to exist already.
    add_dataset_to_entity(
            new_ref_genome, ref_genome_dataset_type, ref_genome_dataset_type,
            new_file_dest)

    return {
        'is_success': True,
        'new_reference_genome': new_ref_genome
    }
