"""
Tests for sv_calling.py
"""

import os

from django.test import TestCase
import vcf

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import get_dataset_with_type
from main.models import Project
from main.models import User
from main.models import Variant
from pipeline.pipeline_runner import run_pipeline_multiple_ref_genomes
from pipeline.snv_calling import get_variant_tool_params
from pipeline.snv_calling import find_variants_with_tool
from scripts.import_util import add_dataset_to_entity
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import copy_dataset_to_entity_data_dir
from scripts.import_util import import_reference_genome_from_local_file
from settings import PWD as GD_ROOT


TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'small_data',
        'ref.fa')

TEST_FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'small_data',
        'simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'small_data',
        'simLibrary.2.fq')

TEST_SAMPLE_UID = '38d786f2'

TEST_BAM = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'small_data',
        'final.bam')

TEST_BAM_INDEX = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'small_data',
        'final.bam.bai')



class TestSVCallers(TestCase):

    def setUp(self):
        user = User.objects.create_user('test_username', password='password',
                email='test@example.com')

        # Grab a project.
        self.project = Project.objects.create(title='test project',
                owner=user.get_profile())

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')

    def test_call_svs(self):
        """Test running the pipeline that finds structural variants.

        The data file consists of 20,000 bases. At 5,000 bases there is
        a 400 base deletion. At 10,000 bases there is a 400 base inversion.
        At 15,000 bases there is a 400 base tandem duplication.

        It seems that Pindel cannot find the inversion. Fortunately,
        delly can usually find inversions; unfortunately, delly only
        works well on large data, so we will not test it here.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Create a sample.
        sample_1 = ExperimentSample.objects.create(
                uid=TEST_SAMPLE_UID,
                project=self.project,
                label='sample1')
        ### Add the raw reads
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)

        # Create relationship between alignment and sample.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=sample_1)
        ### Add alignment data. NOTE: Stored in sample model dir.

        # index (no dataset)
        copy_dataset_to_entity_data_dir(sample_1, TEST_BAM_INDEX)

        # bam file (with dataset)
        copy_dest = copy_dataset_to_entity_data_dir(sample_1, TEST_BAM)
        add_dataset_to_entity(sample_alignment, Dataset.TYPE.BWA_ALIGN,
                Dataset.TYPE.BWA_ALIGN, copy_dest)

        # Make sure there are no variants before.
        self.assertEqual(0, len(Variant.objects.filter(
                reference_genome=self.reference_genome)))

        # Run the pipeline.
        for variant_params in get_variant_tool_params()[1:3]:  # pindel & delly
            find_variants_with_tool(alignment_group, variant_params, project=self.project)

        # Check that the alignment group has a freebayes vcf dataset associated
        # with it.
        vcf_dataset = get_dataset_with_type(alignment_group,
                Dataset.TYPE.VCF_PINDEL)
        self.assertIsNotNone(vcf_dataset)

        # Make sure the .vcf file actually exists.
        self.assertTrue(os.path.exists(vcf_dataset.get_absolute_location()))

        # Make sure the vcf is valid by reading it using pyvcf.
        with open(vcf_dataset.get_absolute_location()) as vcf_fh:
            try:
                reader = vcf.Reader(vcf_fh)
                reader.next()
            except:
                self.fail("Not valid vcf")

        # Grab the resulting variants.
        variants = Variant.objects.filter(reference_genome=self.reference_genome)

        # Check that there is a deletion around base 5000.
        # Check that there is a tandem duplication around base 15000.
        foundDeletion = False
        foundDuplication = False
        print [(variant.position, variant.type) for variant in variants]
        for variant in variants:
            if abs(variant.position - 5000) <= 3 and variant.type == 'DELETION':
                foundDeletion = True
            elif abs(variant.position - 15000) <= 3 and variant.type == 'DUPLICATION':
                foundDuplication = True
        self.assertTrue(foundDeletion)
        self.assertTrue(foundDuplication)


class TestSVPipeline(TestCase):

    def setUp(self):
        user = User.objects.create_user('test_username_sv', password='password',
                email='test@example.com')

        # Grab a project.
        self.project = Project.objects.create(title='test project',
                owner=user.get_profile())

        # Create a ref genome.
        REF = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'all_svs', 'ref.fa')
        FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'all_svs', 'simLibrary.1.fq')
        FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'sv_testing', 'all_svs', 'simLibrary.2.fq')
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', REF, 'fasta')

        self.experiment_sample = ExperimentSample.objects.create(
                project=self.project, label='sample1')
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, FASTQ1)
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, FASTQ2)

    def test_pipeline_and_svs(self):
        run_pipeline_multiple_ref_genomes('name', [self.reference_genome],
                [self.experiment_sample])

        alignment_group_obj_list = AlignmentGroup.objects.filter(
                reference_genome=self.reference_genome)
        self.assertEqual(1, len(alignment_group_obj_list))

        alignment_group_obj = alignment_group_obj_list[0]
        self.assertEqual(1,
                len(alignment_group_obj.experimentsampletoalignment_set.all()))

        # Make sure the initial JBrowse config has been created.
        jbrowse_dir = self.reference_genome.get_jbrowse_directory_path()
        self.assertTrue(os.path.exists(jbrowse_dir))
        self.assertTrue(os.path.exists(os.path.join(jbrowse_dir,
                'indiv_tracks')))

        variant_params = get_variant_tool_params()  # first one is freebayes
        for variant_param in variant_params:
            find_variants_with_tool(alignment_group, variant_param, project=self.project)

        for vcf_type in [Dataset.TYPE.VCF_FREEBAYES,
                Dataset.TYPE.VCF_PINDEL, Dataset.TYPE.VCF_DELLY]:
            print 'checking type', vcf_type
            vcf_dataset = get_dataset_with_type(alignment_group,
                    Dataset.TYPE.VCF_FREEBAYES)
            self.assertIsNotNone(vcf_dataset)
            self.assertTrue(os.path.exists(vcf_dataset.get_absolute_location()))

