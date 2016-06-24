import os
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from main.models import AlignmentGroup
from main.models import Contig


if __name__ == '__main__':
    ag = AlignmentGroup.objects.get(label='Alignment 1')
    esta = ag.experimentsampletoalignment_set.all()[0]

    for c_id in range(300):
        Contig.objects.create(
                label='c_{id}'.format(id=str(c_id)),
                parent_reference_genome=ag.reference_genome,
                experiment_sample_to_alignment=esta)
