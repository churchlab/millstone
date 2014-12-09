"""
Utility methods for testing.
"""

import fileinput
import os
import vcf

from django.contrib.auth.models import User

from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import Dataset

from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSampleToAlignment
from settings import PWD as GD_ROOT
from utils.import_util import copy_and_add_dataset_source
from variants.vcf_parser import parse_alignment_group_vcf


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'
TEST_PROJECT_NAME = 'recoli'
TEST_REF_GENOME_LABEL = 'mg1655'

TEST_GENOME_SNPS = os.path.join(GD_ROOT, 'test_data',
        'fake_genome_and_reads',
        'test_genome_snps.vcf')


def create_common_entities():
    """Creates the most common entities for testing.

    Returns at a User, Project, ReferenceGenome that are all
    related.
    """
    user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

    project = Project.objects.create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())

    reference_genome = ReferenceGenome.objects.create(
            project=project,
            label=TEST_REF_GENOME_LABEL)

    chromosome = Chromosome.objects.create(
            reference_genome=reference_genome,
            label='Chromosome',
            num_bases=9001)

    sample_1 = ExperimentSample.objects.create(
            project=project,
            label='es1')

    sample_2 = ExperimentSample.objects.create(
            project=project,
            label='es2')

    alignment_group_1 = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=reference_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    return {
        'user': user,
        'project': project,
        'reference_genome': reference_genome,
        'chromosome': chromosome,
        'sample_1': sample_1,
        'sample_2': sample_2,
        'alignment_group_1': alignment_group_1
    }

def create_common_entities_w_variants():
    """Creates the most common entities for testing.

    Returns at a User, Project, ReferenceGenome, alignment,
    and variants that are all related.
    """

    # this is the number of samples in the VCF file
    num_samples = 10

    user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

    project = Project.objects.create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())

    reference_genome = ReferenceGenome.objects.create(
            project=project,
            label=TEST_REF_GENOME_LABEL)

    chromosome = Chromosome.objects.create(
            reference_genome=reference_genome,
            label='Chromosome',
            num_bases=9001)

    alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=reference_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    Chromosome.objects.create(
        reference_genome=reference_genome,
        label='Chromosome',
        num_bases=2000)

    VCF_DATATYPE = Dataset.TYPE.VCF_FREEBAYES
    copy_and_add_dataset_source(
                alignment_group, VCF_DATATYPE,
                VCF_DATATYPE, TEST_GENOME_SNPS)

    # Create experiment sample objects having UIDs that correspond to those
    # in the vcf file. This is a bit "fake" in that the actual pipeline we
    # will be generating the vcf file from the samples (see add_groups()
    # stage of pipeline.
    with open(TEST_GENOME_SNPS) as fh:
        reader = vcf.Reader(fh)
        experiment_sample_uids = reader.samples
    samples = [ExperimentSample.objects.create(
                    uid=sample_uid,
                    project=project,
                    label='fakename:' + sample_uid)
            for sample_uid in experiment_sample_uids]

    # add samples to alignment group
    for sample in samples:
        ExperimentSampleToAlignment.objects.get_or_create(
                alignment_group=alignment_group,
                experiment_sample=sample)


    # Parse the vcf
    parse_alignment_group_vcf(alignment_group, VCF_DATATYPE)

    return {
        'user': user,
        'project': project,
        'reference_genome': reference_genome,
        'chromosome': chromosome,
        'samples': samples,
        'alignment_group': alignment_group
    }
