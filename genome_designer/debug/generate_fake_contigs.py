import os
import random
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from main.models import AlignmentGroup
from main.models import Contig
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment


if __name__ == '__main__':
    ag = AlignmentGroup.objects.get(label='Alignment 1')

    esta = ag.experimentsampletoalignment_set.all()[0]

    for c_id in range(300):
        c = Contig.objects.create(
                label='s1_c_{id}'.format(id=str(c_id)),
                parent_reference_genome=ag.reference_genome,
                experiment_sample_to_alignment=esta,
                num_bases=random.randint(0, 100))
        c.metadata['coverage'] = random.random() * 100
        c.save()

    # Create a second sample and set of contigs to be able to test sorting.
    (sample_2, created) = ExperimentSample.objects.get_or_create(
            project=ag.reference_genome.project,
            label='sample2')

    esta_2 = ExperimentSampleToAlignment.objects.create(
            alignment_group=ag,
            experiment_sample=sample_2)

    for c_id in range(300):
        Contig.objects.create(
                label='s2_c_{id}'.format(id=str(c_id)),
                parent_reference_genome=ag.reference_genome,
                experiment_sample_to_alignment=esta_2,
                num_bases=random.randint(0, 100))
        c.metadata['coverage'] = random.random() * 100
        c.save()
