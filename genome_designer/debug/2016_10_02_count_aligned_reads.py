"""
Script to count aligned reads.
"""

import os
import subprocess
import sys

import numpy as np

# Setup Django environment.
sys.path.append(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from django.conf import settings

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.models import Project


def count_alignment_group_reads(ag_uid):
    ag = AlignmentGroup.objects.get(uid=ag_uid)

    esta_set = ag.experimentsampletoalignment_set.all()
    esta_reads = [count_alignment_reads(esta.uid) for esta in esta_set]

    return {
        'num_samples': len(esta_reads),
        'mean_reads': np.mean(esta_reads),
        'stdev_reads': np.std(esta_reads),
    }


def count_alignment_reads(esta_uid):
    esta = ExperimentSampleToAlignment.objects.get(uid=esta_uid)
    bam_path = esta.dataset_set.get(type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()

    cmd = ' '.join([
        settings.SAMTOOLS_BINARY,
        'idxstats',
        bam_path,
        '|',
        'cut', '-f', '1,3'
    ])

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    total_reads = 0
    for line in p.stdout:
        total_reads += int(line.split()[1])
    return total_reads


if __name__ == '__main__':

    for p in Project.objects.all():
        print '>>>>>>>>>>'
        print p.title
        ref_genome = p.referencegenome_set.all()[0]
        ag_set = ref_genome.alignmentgroup_set.all()
        assert len(ag_set) == 1
        ag_set[0].uid
        print count_alignment_group_reads(ag_set[0].uid)
