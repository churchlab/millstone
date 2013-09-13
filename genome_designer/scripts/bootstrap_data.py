#!/usr/bin/env python

"""
Script to setup some test data.

This is useful during development when we are continuously wiping the db
and want to get some new data in quickly.

NOTE: Several tests use this module, so avoid breaking tests when changing
this.
"""

# Since this script is intended to be used from the terminal, setup the
# environment first so that django and model imports work.
from util import setup_django_env
setup_django_env()

import os
import random
import shutil

from django.db import transaction
from django.contrib.auth.models import User
from django.core.management import call_command

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.alignment_pipeline import create_alignment_groups_and_start_alignments
from scripts.snp_callers import run_snp_calling_pipeline
from scripts.import_util import add_dataset_to_entity
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import copy_dataset_to_entity_data_dir
from scripts.import_util import import_reference_genome_from_local_file

import settings
from settings import PWD as GD_ROOT


# This is the directory where this bootstrap script is located.
PWD = os.path.dirname(os.path.realpath(__file__ ))

# Test data.
TEST_USERNAME = 'gmcdev'

TEST_PASSWORD = 'g3n3d3z'

TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

TEST_BAM = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'bwa_align.sorted.grouped.realigned.bam')

TEST_BAM_INDEX = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'bwa_align.sorted.grouped.realigned.bam.bai')

TEST_PROJECT_NAME = 'recoli'

REF_GENOME_1_LABEL = 'mg1655'

REF_GENOME_2_LABEL = 'c321D'

SAMPLE_1_LABEL = 'sample1'

VARIANTSET_1_LABEL = 'Set A'

VARIANTSET_2_LABEL = 'Set B'

# A set of data consisting of a small annotated genome, many samples, and some
# designed SNPs which are each in some of the samples.
class FullVCFTestSet:
    TEST_DIR = os.path.join(GD_ROOT,'test_data', 'full_vcf_test_set')
    NUM_SAMPLES = 2
    TEST_GENBANK = os.path.join(TEST_DIR, 'mg1655_tolC_through_zupT.gb')
    FASTQ1 = [os.path.join(TEST_DIR, 'sample%d.simLibrary.1.fq' % i)
             for i in range(NUM_SAMPLES)]
    FASTQ2 = [os.path.join(TEST_DIR, 'sample%d.simLibrary.2.fq' % i)
             for i in range(NUM_SAMPLES)]
    TEST_DESIGNED_SNPS = os.path.join(TEST_DIR, 'designed_snps.vcf')

def bootstrap_fake_data():
    """Fill the database with fake data.
    """
    ### Get or create the user.
    try:
        user = User.objects.get(username=TEST_USERNAME)
    except User.DoesNotExist:
        user = User.objects.create_user(
                TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

    ### Create some projects
    (test_project, project_created) = Project.objects.get_or_create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())
    (test_project_2, project_created) = Project.objects.get_or_create(
            title='project2', owner=user.get_profile())
    (test_project_3, project_created) = Project.objects.get_or_create(
            title='project3', owner=user.get_profile())

    ### Create some reference genomes
    ref_genome_1 = import_reference_genome_from_local_file(
            test_project, REF_GENOME_1_LABEL, TEST_FASTA, 'fasta')

    ref_genome_2 = import_reference_genome_from_local_file(
            test_project, REF_GENOME_2_LABEL, TEST_FASTA, 'fasta')

    ref_genome_3 = import_reference_genome_from_local_file(
            test_project, 'test_genome', TEST_FASTA, 'fasta')

    ### Create some ExperimentSamples.

    # Create some samples without backing data just to explore the UI.
    ExperimentSample.objects.get_or_create(
            project=test_project,
            label='C321D_MiSeq',
            group='Plate 1',
            well='A01',
            num_reads=25029296,
    )

    ExperimentSample.objects.get_or_create(
            project=test_project,
            label='C321D Fixed 01',
            group='Plate 2',
            well='A01',
            num_reads=12345,
    )

    ExperimentSample.objects.get_or_create(
            project=test_project,
            label='C321D Fixed 02',
            group='Plate 2',
            well='A02',
            num_reads=897213,
    )

    # Create some samples with backing data.
    (sample_1, created) = ExperimentSample.objects.get_or_create(
            project=test_project,
            label=SAMPLE_1_LABEL)
    # Add datasets to the samples.
    if not sample_1.dataset_set.filter(type=Dataset.TYPE.FASTQ1):
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
    if not sample_1.dataset_set.filter(type=Dataset.TYPE.FASTQ2):
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)

    ### Create an alignment.
    alignment_group_1 = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=ref_genome_3,
            aligner=AlignmentGroup.ALIGNER.BWA)
    # Link it to a sample.
    sample_alignment = ExperimentSampleToAlignment.objects.create(
            alignment_group=alignment_group_1,
            experiment_sample=sample_1)
    ### Add alignment data. NOTE: Stored in sample model dir.
    # NOTE: This is a bit convoluted. Perhaps it would be better to store alignments
    # in the ExperimentSampleToAlignment directory.
    copy_dest = copy_dataset_to_entity_data_dir(sample_1, TEST_BAM)
    copy_dataset_to_entity_data_dir(sample_1, TEST_BAM_INDEX)
    add_dataset_to_entity(sample_alignment, Dataset.TYPE.BWA_ALIGN,
            Dataset.TYPE.BWA_ALIGN, copy_dest)

    ### Create some fake variants
    @transaction.commit_on_success
    def _create_fake_variants():
        for var_count in range(100):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=ref_genome_1,
                chromosome='chrom',
                position=random.randint(1,ref_genome_1.num_bases),
                ref_value='A',
                alt_value='G')
    _create_fake_variants()

    ### Add fake variants to a set
    @transaction.commit_on_success
    def _add_fake_variants_to_fake_set():
        ref_genome_1 = ReferenceGenome.objects.get(
            label=REF_GENOME_1_LABEL)

        (sample_1, created) = ExperimentSample.objects.get_or_create(
            project=test_project,
            label=SAMPLE_1_LABEL)

        var_set1 = VariantSet.objects.create(
            reference_genome=ref_genome_1,
            label=VARIANTSET_1_LABEL)
        var_set2 = VariantSet.objects.create(
            reference_genome=ref_genome_1,
            label=VARIANTSET_2_LABEL)

        variant_list = Variant.objects.filter(
            reference_genome=ref_genome_1)
        for var in variant_list:

            #add variant to one of two sets, depending on var position
            if var.position < 50:
                if var.position < 25:
                    vvs1 = VariantToVariantSet.objects.create(
                        variant=var,
                        variant_set=var_set1)

                    #add a sample to the association if the variant is odd
                    if var.position % 2:
                        vvs1.sample_variant_set_association.add(sample_1)

                if var.position > 20:
                    vvs2 = VariantToVariantSet.objects.create(
                        variant=var,
                        variant_set=var_set2)

                    #add a sample to the association if the variant is even
                    if not var.position % 2:
                        vvs2.sample_variant_set_association.add(sample_1)

        # Make sure that both sets have at least one variant.
        guaranteed_var = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=ref_genome_1,
                chromosome='chrom',
                position=22,
                ref_value='A',
                alt_value='G')
        vvs1 = VariantToVariantSet.objects.create(
                variant=guaranteed_var,
                variant_set=var_set1)
        vvs2 = VariantToVariantSet.objects.create(
                variant=guaranteed_var,
                variant_set=var_set2)
        vvs2.sample_variant_set_association.add(sample_1)

    _add_fake_variants_to_fake_set()

    #############################
    # Full VCF Testing (annotated for snpeff, variant filtering, etc)
    #############################

    # Create a new reference genome and samples using full_vcf_test_set
    full_vcf_reference_genome = import_reference_genome_from_local_file(
                test_project, 'mg1655_tolC_through_zupT',
                FullVCFTestSet.TEST_GENBANK, 'genbank')

    # Create all samples.
    full_vcf_samples = []
    for i in range(FullVCFTestSet.NUM_SAMPLES):
        sample_obj = ExperimentSample.objects.create(
                project=test_project,
                label='Sample %d' % i)

        # Add raw reads to each sample.
        copy_and_add_dataset_source(sample_obj, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, FullVCFTestSet.FASTQ1[i])
        copy_and_add_dataset_source(sample_obj, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, FullVCFTestSet.FASTQ2[i])

        full_vcf_samples.append(sample_obj)

    # Run the alignment. Return the alignment group created, indexed by the
    #   reference genome's uid.
    full_vcf_alignment_group = create_alignment_groups_and_start_alignments(
           [full_vcf_reference_genome],
           full_vcf_samples)[full_vcf_reference_genome.uid]


    # Run Freebayes. This should also kick off snpeff afterwards.
    run_snp_calling_pipeline(full_vcf_alignment_group)



def reset_database():
    """Deletes the old database and sets up a new one.

    For now, only works with the temp.db database to prevent
    accidentally deleting data down the line.
    """
    ### Delete the old database if it exists.
    print 'Deleting old database ...'
    TEMP_DB_NAME = 'temp.db'
    temp_db_name = settings.DATABASES['default']['NAME']
    assert temp_db_name == TEMP_DB_NAME
    tempdb_dir = os.path.split(PWD)[0]
    tempdb_abs_path = os.path.join(tempdb_dir, TEMP_DB_NAME)
    if os.path.exists(tempdb_abs_path):
        os.remove(tempdb_abs_path)

    ### Run syncdb
    # NOTE: Remove interactive=False if you want to have the option of creating
    # a super user on sync.
    print 'Creating new database via syncdb ...'
    call_command('syncdb', interactive=False)

    ### Recreate the media root.
    if os.path.exists(settings.MEDIA_ROOT):
        shutil.rmtree(settings.MEDIA_ROOT)
    os.mkdir(settings.MEDIA_ROOT)


if __name__ == '__main__':
    reset_database()
    bootstrap_fake_data()
