#!/usr/bin/env python

"""
Script to setup some test data.

This is useful during development when we are continuously wiping the db
and want to get some new data in quickly.

NOTE: Several tests use this module, so avoid breaking tests when changing
this.
"""

import os
import random
import shutil
import subprocess
import sys

# Setup the Django environment. Our setup_django_env() method is outdated.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from django.db import transaction
from django.contrib.auth.models import User
from django.core.management import call_command
from django.conf import settings

from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import Region
from main.models import RegionInterval
from main.models import SavedVariantFilterQuery
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import FullVCFTestSet
from pipeline.pipeline_runner import run_pipeline
from utils.import_util import add_dataset_to_entity
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import copy_dataset_to_entity_data_dir
from utils.import_util import import_reference_genome_from_local_file
from utils.import_util import import_variant_set_from_vcf
from utils.import_util import run_fastqc_on_sample_fastq

from settings import PWD as GD_ROOT

# Test data.
TEST_USERNAME = 'gmcdev'

TEST_PASSWORD = 'g3n3d3z'

TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_FASTA = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

TEST_FASTQ_GZ_1 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '6057f443', 'test_genome_8.snps.simLibrary.1.fq.gz')

TEST_FASTQ_GZ_2 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '6057f443', 'test_genome_8.snps.simLibrary.2.fq.gz')

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

CUSTOM_SAVED_QUERY_LIST = [
    'GT_TYPE = 2 & DP > 10',
    'INFO_EFF_IMPACT = HIGH',
]

@transaction.commit_on_success
def add_fake_variants_to_genome(ref_genome, chromosome, num_variants):
    """
    Add many fake variants to a fake genome in order to test
    features like pagination.
    """

    for var_count in range(num_variants):

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=ref_genome,
                chromosome=chromosome,
                position=random.randint(1, chromosome.num_bases),
                ref_value='A')

        variant.variantalternate_set.add(VariantAlternate.objects.create(
                variant=variant,
                alt_value='G'))

def create_fake_variants_and_variant_sets(ref_genome, num_variants= 100):

    chrom_0 = Chromosome.objects.get(reference_genome=ref_genome)

    add_fake_variants_to_genome(
            ref_genome=ref_genome,
            chromosome=chrom_0,
            num_variants=num_variants)

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
                chromosome=Chromosome.objects.get(
                    reference_genome=ref_genome),
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
        print "trying to create user"
        return User.objects.get(username=TEST_USERNAME)
    except User.DoesNotExist:
        print "user did not exist"
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

    ### Create some saved queries.
    for saved_query_text in CUSTOM_SAVED_QUERY_LIST:
        SavedVariantFilterQuery.objects.get_or_create(
                owner=user.get_profile(),
                text=saved_query_text)

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

    # Create sample backed by g-zipped data.
    gz_backed_sample = ExperimentSample.objects.create(
            project=test_project,
            label='sample backed by gz data')
    gz_fastq1_dataset = copy_and_add_dataset_source(
            gz_backed_sample, Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ1,
            TEST_FASTQ_GZ_1)
    gz_fastq2_dataset = copy_and_add_dataset_source(
            gz_backed_sample, Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ2,
            TEST_FASTQ_GZ_2)
    run_fastqc_on_sample_fastq(gz_backed_sample, gz_fastq1_dataset)
    run_fastqc_on_sample_fastq(gz_backed_sample, gz_fastq2_dataset, rev=True)

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
    parent_obj = None
    full_vcf_samples = []
    for i in range(FullVCFTestSet.NUM_SAMPLES):
        sample_obj = ExperimentSample.objects.create(
                project=test_project,
                label='Sample %d' % i)

        sample_obj.data['SAMPLE_WELL'] = 'A0%d' % (i+1)

        if i == 0:
            parent_obj = sample_obj
        if i > 0:
            sample_obj.data['SAMPLE_PARENTS'] = parent_obj.label
            parent_obj.add_child(sample_obj)
            parent_obj.save()

        sample_obj.save()

        # Add raw reads to each sample.
        fastq1_dataset = copy_and_add_dataset_source(sample_obj,
                Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1,
                FullVCFTestSet.FASTQ1[i])
        fastq2_dataset = copy_and_add_dataset_source(sample_obj,
                Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2,
                FullVCFTestSet.FASTQ2[i])

        # Run FASTQC on sample reads.
        run_fastqc_on_sample_fastq(sample_obj, fastq1_dataset)
        run_fastqc_on_sample_fastq(sample_obj, fastq2_dataset, rev=True)

        full_vcf_samples.append(sample_obj)

    # Run the alignment. Return the alignment group created, indexed by the
    # reference genome's uid.
    run_pipeline(
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

    script_string = """
    sudo -u %(os_user)s psql -c "
    DO
    \$body\$
    BEGIN
       IF NOT EXISTS (
          SELECT *
          FROM   pg_catalog.pg_user
          WHERE  usename = '%(user)s') THEN
          CREATE USER %(user)s WITH PASSWORD '%(password)s';
       END IF;
    END;
    \$body\$"
    sudo -u %(os_user)s psql -c "DROP DATABASE IF EXISTS %(db)s;"
    sudo -u %(os_user)s psql -c "DROP USER IF EXISTS %(user)s;"
    sudo -u %(os_user)s psql -c "CREATE USER %(user)s WITH PASSWORD \
            '%(password)s';"
    sudo -u %(os_user)s psql -c "CREATE DATABASE %(db)s;"
    sudo -u %(os_user)s psql -c 'GRANT ALL PRIVILEGES ON DATABASE \
            %(db)s to %(user)s;'
    sudo -u %(os_user)s psql -c "ALTER USER %(user)s CREATEDB;"
    """ % {
            'db': settings.DATABASES['default']['NAME'],
            'user': settings.DATABASES['default']['USER'],
            'password': settings.DATABASES['default']['PASSWORD'],
            'os_user': settings.DATABASES['default']['OS_USER']
    }

    proc = subprocess.Popen(script_string, shell=True, stderr=subprocess.PIPE)

    # parse script stderr for errors.
    error_lines = []
    for output_line in proc.stderr.readline():
        # Skip 'NOTICE' stderr lines, they are not errors.
        if not output_line or 'NOTICE' in output_line:
            continue
        error_lines.append(output_line)

    if len(error_lines):
        raise Exception(
                'Error while reseting database. Possible reasons:\n '
                '\t* Celery is running\n'
                '\t* Postgres session is open\n'
                '\t* No postgres user; change OS_USER in default database'
                        ' in local_settings.py\n'
                '\nOffending postgres errors:\n' + ''.join(error_lines))

    """
    flush: removes all rows in the database.
    syncdb --all: performs syncdb on non-South apps and migrate on South apps
    migrate --fake: note that migrate has been performed in syncdb
    """
    call_command('syncdb', migrate_all=True, interactive=False)
    print "\n\n\n\n---------------ABOUT TO MIGRATE"
    call_command('migrate', fake=True, interactive=False)
    print "\n\n\n\n----------------MIGRATION FINISHED"

    ### Recreate the media root.
    if os.path.exists(settings.MEDIA_ROOT):
        shutil.rmtree(settings.MEDIA_ROOT)
    os.mkdir(settings.MEDIA_ROOT)


def confirm_bootstrap():
    if len(sys.argv) > 1 and sys.argv[1] in ["-q",]:
        sys.argv.pop(1)
        return True
    confirm_text = raw_input(
            "This will wipe any current database. Are you sure? y/n\n")
    return confirm_text.lower() in ['y', 'yes']


if __name__ == '__main__':

    # First entry is the default.
    # We can add more here later.
    BOOTSTRAP_FLAVORS = [
        ('blank', lambda: 'no function'),
        ('full', bootstrap_fake_data)]

    if confirm_bootstrap():
        if len(sys.argv) == 1:

            reset_database()
            BOOTSTRAP_FLAVORS[0][1]()

        elif len(sys.argv) >= 2:

            # Make sure first arg is an available flavor.
            assert sys.argv[1] in dict(BOOTSTRAP_FLAVORS), (
                    'Available flavors:\n\t%s') % ('\n\t'.join(
                    dict(BOOTSTRAP_FLAVORS).keys()))

            # Reset the database and bootstrap with the flavor.
            reset_database()
            dict(BOOTSTRAP_FLAVORS)[sys.argv[1]]()

    else:
        print 'Aborting.'
