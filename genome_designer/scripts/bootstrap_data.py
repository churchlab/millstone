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
import sys

from django.db import transaction
from django.contrib.auth.models import User
from django.core.management import call_command
from django.conf import settings

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import Region
from main.models import RegionInterval
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from pipeline.pipeline_runner import run_pipeline
from scripts.import_util import add_dataset_to_entity
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import copy_dataset_to_entity_data_dir
from scripts.import_util import import_reference_genome_from_local_file
from scripts.import_util import import_variant_set_from_vcf

from settings import PWD as GD_ROOT

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

SV_PROJECT_NAME = 'sv_testing'

REF_GENOME_1_LABEL = 'mg1655'

REF_GENOME_2_LABEL = 'c321D'

SAMPLE_1_LABEL = 'sample1'

VARIANTSET_1_LABEL = 'Set A'

VARIANTSET_2_LABEL = 'Set B'

# A set of data consisting of a small annotated genome, many samples, and some
# designed SNPs which are each in some of the samples.
class FullVCFTestSet:
    TEST_DIR = os.path.join(GD_ROOT,'test_data', 'full_vcf_test_set')
    NUM_SAMPLES = 3
    TEST_GENBANK = os.path.join(TEST_DIR, 'mg1655_tolC_through_zupT.gb')
    FASTQ1 = [os.path.join(TEST_DIR, 'sample%d.simLibrary.1.fq' % i)
             for i in range(NUM_SAMPLES)]
    FASTQ2 = [os.path.join(TEST_DIR, 'sample%d.simLibrary.2.fq' % i)
             for i in range(NUM_SAMPLES)]
    TEST_DESIGNED_SNPS = os.path.join(TEST_DIR, 'designed_snps.vcf')


def create_fake_variants_and_variant_sets(ref_genome):
    @transaction.commit_on_success
    def _create_fake_variants():
        for var_count in range(100):
            variant = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=ref_genome,
                    chromosome='chrom',
                    position=random.randint(1,ref_genome.num_bases),
                    ref_value='A')

            variant.variantalternate_set.add(VariantAlternate.objects.create(
                    variant=variant,
                    alt_value='G'))

    _create_fake_variants()

    ### Add fake variants to a set
    @transaction.commit_on_success
    def _add_fake_variants_to_fake_set():
        ref_genome_1 = ReferenceGenome.objects.get(
            label=REF_GENOME_1_LABEL)

        (sample_1, created) = ExperimentSample.objects.get_or_create(
            project=ref_genome.project,
            label=SAMPLE_1_LABEL)

        var_set1 = VariantSet.objects.create(
            reference_genome=ref_genome,
            label=VARIANTSET_1_LABEL)
        var_set2 = VariantSet.objects.create(
            reference_genome=ref_genome,
            label=VARIANTSET_2_LABEL)

        variant_list = Variant.objects.filter(reference_genome=ref_genome)
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
                reference_genome=ref_genome,
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


def get_or_create_user():
    try:
        return User.objects.get(username=TEST_USERNAME)
    except User.DoesNotExist:
        return User.objects.create_user(
                TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)


def bootstrap_fake_data():
    """Fill the database with fake data.
    """
    user = get_or_create_user()

    ### Create some projects
    (test_project, project_created) = Project.objects.get_or_create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())
    (test_project_2, project_created) = Project.objects.get_or_create(
            title=SV_PROJECT_NAME, owner=user.get_profile())

    ### Create some reference genomes
    ref_genome_1 = import_reference_genome_from_local_file(
            test_project, REF_GENOME_1_LABEL, TEST_FASTA, 'fasta')

    ref_genome_2 = import_reference_genome_from_local_file(
            test_project, REF_GENOME_2_LABEL, TEST_FASTA, 'fasta')

    ref_genome_3 = import_reference_genome_from_local_file(
            test_project, 'test_genome', TEST_FASTA, 'fasta')

    ### Create some ExperimentSamples.

    # Create some samples without backing data just to explore the UI.
    ExperimentSample.objects.create(
            project=test_project,
            label='C321D_MiSeq',
            data = {'SAMPLE_WELL': 'A01'}
    )

    ExperimentSample.objects.create(
            project=test_project,
            label='C321D Fixed 01',
            data = {'SAMPLE_WELL': 'A02'}
    )

    ExperimentSample.objects.create(
            project=test_project,
            label='C321D Fixed 02',
            data = {'SAMPLE_WELL': 'A03'}
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

    # Create fake variants.
    create_fake_variants_and_variant_sets(ref_genome_1)

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

        sample_obj.data['SAMPLE_WELL'] = 'A0%d' % (i+1)
        sample_obj.save()

        # Add raw reads to each sample.
        copy_and_add_dataset_source(sample_obj, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, FullVCFTestSet.FASTQ1[i])
        copy_and_add_dataset_source(sample_obj, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, FullVCFTestSet.FASTQ2[i])

        full_vcf_samples.append(sample_obj)

    # Run the alignment. Return the alignment group created, indexed by the
    # reference genome's uid.
    full_vcf_alignment_group = run_pipeline(
            'test_align', full_vcf_reference_genome, full_vcf_samples)

    import_variant_set_from_vcf(full_vcf_reference_genome, 'Designed',
            FullVCFTestSet.TEST_DESIGNED_SNPS)

    def _create_region_intervals(region, interval_tuple_list):
        """Helper method to create RegionIntervals for a Region.

        Args:
            region: Region Model object.
            interval_tuple_list: List of tuples of intervals to create.
        """
        for interval in interval_tuple_list:
            RegionInterval.objects.create(
                    region=region,
                    start=interval[0],
                    end=interval[1])

    # Create some fake regions.
    # TODO: Should not be much harder to replace this with real regions.
    region_1 = Region.objects.create(
        reference_genome=full_vcf_reference_genome,
        label='region_1',
        type=Region.TYPE.POOR_MAPPING_QUALITY)
    _create_region_intervals(region_1, [(1,150), (300, 400), (500, 900)])

    region_2 = Region.objects.create(
        reference_genome=full_vcf_reference_genome,
        label='region_2',
        type=Region.TYPE.POOR_MAPPING_QUALITY)
    _create_region_intervals(region_2, [(1000, 1500)])

    region_3 = Region.objects.create(
        reference_genome=full_vcf_reference_genome,
        label='region_3',
        type=Region.TYPE.POOR_MAPPING_QUALITY)
    _create_region_intervals(region_3, [(1800, 1900), (2150, 2300)])

    # And some GENE regions.

    gene_A = Region.objects.create(
        reference_genome=full_vcf_reference_genome,
        label='geneA',
        type=Region.TYPE.GENE)
    _create_region_intervals(gene_A, [(2000, 2400)])

    gene_B = Region.objects.create(
        reference_genome=full_vcf_reference_genome,
        label='geneB',
        type=Region.TYPE.GENE)
    _create_region_intervals(gene_B, [(4800, 5200)])

    gene_C = Region.objects.create(
        reference_genome=full_vcf_reference_genome,
        label='geneC',
        type=Region.TYPE.GENE)
    _create_region_intervals(gene_C, [(1, 500)])

    # Bootstrap test_project_2 with SV stuff
    sv_testing_bootstrap(test_project_2)


def sv_testing_bootstrap(project):
    sv_testing_dir = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'all_svs')
    fasta = os.path.join(sv_testing_dir, 'ref.fa')
    fq1 = os.path.join(sv_testing_dir, 'simLibrary.1.fq')
    fq2 = os.path.join(sv_testing_dir, 'simLibrary.2.fq')

    ref_genome = import_reference_genome_from_local_file(
            project, 'ref', fasta, 'fasta')

    sample = ExperimentSample.objects.create(
            project=project,
            label='simLibrary',
    )
    copy_and_add_dataset_source(sample, Dataset.TYPE.FASTQ1,
            Dataset.TYPE.FASTQ1, fq1)
    copy_and_add_dataset_source(sample, Dataset.TYPE.FASTQ2,
            Dataset.TYPE.FASTQ2, fq2)

    if '--sv' in sys.argv:  # using --sv argument runs pipeline for SV project
        run_pipeline('sample_alignment_ref', ref_genome, [sample])


def reset_database():
    """Deletes the old database and sets up a new one.

    For now, only works with the temp.db database to prevent
    accidentally deleting data down the line.
    """
    ### Delete the old database if it exists.
    print 'Deleting old database ...'

    """
    flush: removes all rows in the database.
    syncdb --all: performs syncdb on non-South apps and migrate on South apps
    migrate --fake: note that migrate has been performed in syncdb
    """
    call_command('flush', interactive=False)
    call_command('syncdb', migrate_all=True, interactive=False)
    call_command('migrate', fake=True, interactive=False)

    ### Recreate the media root.
    if os.path.exists(settings.MEDIA_ROOT):
        shutil.rmtree(settings.MEDIA_ROOT)
    os.mkdir(settings.MEDIA_ROOT)


def confirm_bootstrap():
    if len(sys.argv) > 1 and sys.argv[1] in ["-q",]:
        return True
    confirm_text = raw_input(
            "This will wipe any current database. Are you sure? y/n\n")
    return confirm_text.lower() in ['y', 'yes']


if __name__ == '__main__':
    if confirm_bootstrap():
        reset_database()
        bootstrap_fake_data()
    else:
        print 'Aborting'
