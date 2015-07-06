"""
These tests cover the individual code pieces related to structural variant
calling.
"""

import os

from django.conf import settings
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
from main.models import VariantAlternate
from pipeline.variant_calling import find_variants_with_tool
from pipeline.variant_calling import VARIANT_TOOL_PARAMS_MAP
from utils.import_util import add_dataset_to_entity
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import copy_dataset_to_entity_data_dir
from utils.import_util import import_reference_genome_from_local_file


TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'sv_testing', 'small_data',
        'ref.fa')

TEST_FASTQ1 = os.path.join(settings.PWD, 'test_data', 'sv_testing', 'small_data',
        'simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(settings.PWD, 'test_data', 'sv_testing', 'small_data',
        'simLibrary.2.fq')

TEST_SAMPLE_UID = '38d786f2'

TEST_BAM = os.path.join(settings.PWD, 'test_data', 'sv_testing', 'small_data',
        'final.bam')

TEST_BAM_INDEX = os.path.join(settings.PWD, 'test_data', 'sv_testing', 'small_data',
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

    def test_end_to_end(self):
        """Test running full pipline on small-ish data.

        The data file consists of 20,000 bases. At 5,000 bases there is
        a 400 base deletion. At 10,000 bases there is a 400 base inversion.
        At 15,000 bases there is a 400 base tandem duplication.

        It seems that Pindel cannot find the inversion. Fortunately,
        delly can usually find inversions. Unfortunately, delly only
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

        # Test with Pindel only for now.
        for tool in ['pindel']:
            find_variants_with_tool(alignment_group,
                    VARIANT_TOOL_PARAMS_MAP[tool], project=self.project)

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

        # Confirm that 2 variants found.
        self.assertEqual(2, len(variants))

        variant_map = {}
        for variant in variants:
            variant_alternates = VariantAlternate.objects.filter(variant=variant)

            # There should be only one variant alternate per SV.
            self.assertEqual(len(variant_alternates), 1)

            pos = variant.position
            svtype = variant_alternates[0].data['INFO_SVTYPE']
            svlen = variant_alternates[0].data['INFO_SVLEN']
            variant_map[svtype] = (pos, svlen)

        # Check that there is a deletion around base 5000.
        self.assertTrue('DEL' in variant_map)
        self.assertTrue(abs(variant_map['DEL'][0] - 5000) <= 3)
        self.assertTrue(abs(variant_map['DEL'][1] - 400) <= 3)

        # Check that there is a tandem duplication around base 15000.
        self.assertTrue('DUP:TANDEM' in variant_map)
        self.assertTrue(abs(variant_map['DUP:TANDEM'][0] - 15000) <= 3)
        self.assertTrue(abs(variant_map['DUP:TANDEM'][1] - 400) <= 3)
