import os
import random
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from main.models import AlignmentGroup
from main.models import Contig
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from utils.import_util import add_dataset_to_entity


def _make_fake_contig(label, esta):
    c = Contig.objects.create(
            label=label,
            parent_reference_genome=ag.reference_genome,
            experiment_sample_to_alignment=esta)
    c.metadata['coverage'] = random.random() * 100

    # Add fasta.
    c.ensure_model_data_dir_exists()

    # Random sequence.
    num_bases = random.randint(0, 100)
    seq = Seq(''.join([random.choice('ATCG') for i in range(num_bases)]))
    seq_record = SeqRecord(seq, id=c.uid)
    dataset_path = os.path.join(c.get_model_data_dir(), 'fasta.fa')
    with open(dataset_path, 'w') as fh:
        SeqIO.write(seq_record, fh, 'fasta')
    add_dataset_to_entity(
            c,
            'contig_fasta',
            Dataset.TYPE.REFERENCE_GENOME_FASTA,
            filesystem_location=dataset_path)

    c.save()


if __name__ == '__main__':
    ag = AlignmentGroup.objects.get(label='Alignment 1')

    esta = ag.experimentsampletoalignment_set.all()[0]

    for c_id in range(300):
        _make_fake_contig('s1_c_{id}'.format(id=str(c_id)), esta)

    # Create a second sample and set of contigs to be able to test sorting.
    (sample_2, created) = ExperimentSample.objects.get_or_create(
            project=ag.reference_genome.project,
            label='sample2')

    esta_2 = ExperimentSampleToAlignment.objects.create(
            alignment_group=ag,
            experiment_sample=sample_2)

    for c_id in range(300):
        _make_fake_contig('s2_c_{id}'.format(id=str(c_id)), esta_2)
