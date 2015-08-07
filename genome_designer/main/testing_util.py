"""
Utility methods for testing.
"""

import os
import vcf

from Bio import SeqIO
from django.conf import settings
from django.contrib.auth.models import User

from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSampleToAlignment
from utils.import_util import copy_and_add_dataset_source
from variants.vcf_parser import parse_alignment_group_vcf
from variants.vcf_parser import parse_vcf


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'
TEST_PROJECT_NAME = 'recoli'
TEST_REF_GENOME_LABEL = 'mg1655'

TEST_DATA_DIR = os.path.join(settings.PWD, 'test_data')

TEST_GENOME_SNPS = os.path.join(settings.PWD, 'test_data',
        'fake_genome_and_reads',
        'test_genome_snps.vcf')


# A set of data consisting of a small annotated genome, many samples, and some
# designed SNPs which are each in some of the samples.
class FullVCFTestSet:
    TEST_DIR = os.path.join(settings.PWD, 'test_data', 'full_vcf_test_set')
    NUM_SAMPLES = 3
    TEST_GENBANK = os.path.join(TEST_DIR, 'mg1655_tolC_through_zupT.gb')
    FASTQ1 = [os.path.join(TEST_DIR, 'sample%d.simLibrary.1.fq' % i)
             for i in range(NUM_SAMPLES)]
    FASTQ2 = [os.path.join(TEST_DIR, 'sample%d.simLibrary.2.fq' % i)
             for i in range(NUM_SAMPLES)]
    TEST_DESIGNED_SNPS = os.path.join(TEST_DIR, 'designed_snps.vcf')


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
            seqrecord_id='Chromosome',
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


def create_sample_and_alignment(
        project, alignment_group, sample_uid, bwa_alignment=None):
    sample = ExperimentSample.objects.create(
            uid=sample_uid, project=project, label=sample_uid)
    sample_alignment = ExperimentSampleToAlignment.objects.create(
            alignment_group=alignment_group, experiment_sample=sample)
    if bwa_alignment is not None:
        copy_and_add_dataset_source(
                sample_alignment, Dataset.TYPE.BWA_ALIGN, Dataset.TYPE.BWA_ALIGN,
                bwa_alignment)
    return {
        'sample': sample,
        'sample_alignment': sample_alignment
    }


def create_recoli_sv_data_from_vcf(project):
    """Populate database with SVs from lumpy vcf.
    """
    VCF_PARSER_TEST_DATA_DIR = os.path.join(TEST_DATA_DIR, 'vcf_parser_test_data')

    LUMPY_4_SAMPLES_RECOLI_VCF = os.path.join(
            VCF_PARSER_TEST_DATA_DIR, 'lumpy_4_samples_recoli.vcf')

    SAMPLE_1_UID = '3990b0f4'
    SAMPLE_2_UID = '0e250e34'
    SAMPLE_3_UID = '396ea926'
    SAMPLE_4_UID = '4a09d3dd'

    reference_genome = ReferenceGenome.objects.create(
            project=project, label='myref')

    Chromosome.objects.create(
            reference_genome=reference_genome,
            label='the chrom',
            seqrecord_id='U00096.2',
            num_bases=5000000000)

    alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1', reference_genome=reference_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    # Connect lumpy vcf as Dataset.
    lumpy_vcf_dataset = copy_and_add_dataset_source(
            alignment_group, Dataset.TYPE.VCF_LUMPY, Dataset.TYPE.VCF_LUMPY,
            LUMPY_4_SAMPLES_RECOLI_VCF)

    # Create samples corresponding to sample ids in vcf.
    create_sample_and_alignment(
            project, alignment_group, SAMPLE_1_UID)
    create_sample_and_alignment(
            project, alignment_group, SAMPLE_2_UID)
    create_sample_and_alignment(
            project, alignment_group, SAMPLE_3_UID)
    create_sample_and_alignment(
            project, alignment_group, SAMPLE_4_UID)

    # Now we have everything we need to parse the vcf.
    parse_vcf(lumpy_vcf_dataset, alignment_group)


def are_fastas_same(fasta_1, fasta_2):
    """"Returns tuple ((bool)fastas equal, (list)indexes of dissimilarity)
    """

    with open(fasta_1, 'r') as fasta_1_fh, \
         open(fasta_2, 'r') as fasta_2_fh:

        fasta_1_seqrecord_list = list(SeqIO.parse(fasta_1_fh, 'fasta'))
        fasta_2_seqrecord_list = list(SeqIO.parse(fasta_2_fh, 'fasta'))

        assert len(fasta_1_seqrecord_list) == 1
        assert len(fasta_2_seqrecord_list) == 1

        seq_1 = fasta_1_seqrecord_list[0].seq
        seq_2 = fasta_2_seqrecord_list[0].seq

        if str(seq_1) == str(seq_2):
            return (True, [])
        else:
            eq = map(lambda x, y: x == y, seq_1, seq_2)
            indexes = [i for i, x in enumerate(eq) if x == 0]
            return (False, indexes)
