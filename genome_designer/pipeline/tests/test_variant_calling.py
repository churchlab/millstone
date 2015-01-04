"""
Tests for variant_calling.py
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
from main.model_utils import clean_filesystem_location
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_split_reads
from pipeline.variant_calling import find_variants_with_tool
from pipeline.variant_calling import run_lumpy
from pipeline.variant_calling import VARIANT_TOOL_PARAMS_MAP
from pipeline.variant_calling.freebayes import freebayes_regions
from pipeline.variant_calling.freebayes import merge_freebayes_parallel
from settings import PWD as GD_ROOT
from utils.import_util import add_dataset_to_entity
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import copy_dataset_to_entity_data_dir
from utils.import_util import import_reference_genome_from_local_file
from variants.vcf_parser import parse_alignment_group_vcf


TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

TEST_SAMPLE_UID = '38d786f2'

TEST_BAM = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'bwa_align.sorted.grouped.realigned.bam')

TEST_BAM_INDEX = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'bwa_align.sorted.grouped.realigned.bam.bai')

TEST_DISC_SPLIT_BAM = os.path.join(GD_ROOT, 'test_data',
        'discordant_split_reads',
        'bwa_align.bam')

LUMPY_VCF = os.path.join(GD_ROOT, 'test_data',
        'discordant_split_reads',
        'lumpy.vcf')


class TestSNPCallers(TestCase):

    def setUp(self):
        user = User.objects.create_user('test_username', password='password',
                email='test@example.com')

        # Grab a project.
        self.project = Project.objects.create(title='test project',
                owner=user.get_profile())

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')

        self.EXPECTED_NUM_VARIANTS = 30

        # We created the test genome with these specs.
        self.EXPECTED_VARIANT_POSITIONS = [800] # one-indexed.
        while len(self.EXPECTED_VARIANT_POSITIONS) < self.EXPECTED_NUM_VARIANTS:
            self.EXPECTED_VARIANT_POSITIONS.append(
                    self.EXPECTED_VARIANT_POSITIONS[-1] + 20)

        self.KNOWN_SUBSTITUTIONS_ROOT = os.path.join(GD_ROOT, 'test_data',
                'test_genome_known_substitutions')

        self.TEST_GENOME_FASTA = os.path.join(self.KNOWN_SUBSTITUTIONS_ROOT,
                'test_genome_known_substitutions.fa')

        self.FAKE_READS_FASTQ1 = os.path.join(self.KNOWN_SUBSTITUTIONS_ROOT,
                'test_genome_known_substitutions_0.snps.simLibrary.1.fq')

        self.FAKE_READS_FASTQ2 = os.path.join(self.KNOWN_SUBSTITUTIONS_ROOT,
                'test_genome_known_substitutions_0.snps.simLibrary.2.fq')

        self.FAKE_READS_SAMPLE_UID = '93b68da4'

        self.FAKE_READS_BAM = os.path.join(self.KNOWN_SUBSTITUTIONS_ROOT,
                'bwa_align.sorted.grouped.realigned.bam')

        self.FAKE_READS_BAM_INDEX = os.path.join(self.KNOWN_SUBSTITUTIONS_ROOT,
                'bwa_align.sorted.grouped.realigned.bam.bai')

        # Create a ref genome from the above.
        self.REFERENCE_GENOME = import_reference_genome_from_local_file(
                self.project, 'test_genome', self.TEST_GENOME_FASTA, 'fasta')

    def _create_alignment(self, haploid=False):

        # Create a new alignment group.
        self.alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.REFERENCE_GENOME)

        if haploid:
            self.alignment_group.alignment_options['call_as_haploid'] = True

        # Create a sample.
        self.sample_1 = ExperimentSample.objects.create(
                uid=self.FAKE_READS_SAMPLE_UID,
                project=self.project,
                label='sample1')
        ### Add the raw reads
        copy_and_add_dataset_source(self.sample_1, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, self.FAKE_READS_FASTQ1)
        copy_and_add_dataset_source(self.sample_1, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, self.FAKE_READS_FASTQ2)

        # Create alignment to the sample.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=self.alignment_group,
                experiment_sample=self.sample_1)
        ### Add alignment data. NOTE: Stored in sample model dir.
        copy_dest = copy_dataset_to_entity_data_dir(
                self.sample_1, self.FAKE_READS_BAM)
        copy_dataset_to_entity_data_dir(self.sample_1,
                self.FAKE_READS_BAM_INDEX)
        add_dataset_to_entity(sample_alignment, Dataset.TYPE.BWA_ALIGN,
                Dataset.TYPE.BWA_ALIGN, copy_dest)

    def test_freebayes_regions(self):
        regions = freebayes_regions(self.reference_genome, region_size=100)
        self.assertEqual(len(regions), 20)

    def test_call_snvs(self):
        """Test running the pipeline that calls SNPS.

        This test doesn't check the accuracy of the SNP-calling. The test is
        intended just to run the pipeline and make sure there are no errors.
        """

        # create the alignment
        self._create_alignment(haploid=False)

        # Run the pipeline.
        variant_params = VARIANT_TOOL_PARAMS_MAP['freebayes']
        find_variants_with_tool(self.alignment_group, variant_params,
                project=self.project)

        # Check that the alignment group has a freebayes vcf dataset associated
        # with it.
        vcf_dataset = get_dataset_with_type(self.alignment_group,
                Dataset.TYPE.VCF_FREEBAYES)
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


    def _freebayes_tester(self, haploid=False):
        """
        Test that freebayes with the default settings works for control
        input data.

        The test data is generated by starting with a short reference genome,
        manually introducing SNPs to create an alternate reference, and
        then running simNGS to create fake sequencing reads. Finally we align
        the fake sequencing reads back to the unmodified genome.  We expect
        freebayes to nail all of these "obvious" cases.
        """

        # Make sure there are no variants before.
        self.assertEqual(0, len(Variant.objects.filter(
                reference_genome=self.REFERENCE_GENOME)))

        # Run the pipeline.
        variant_params = VARIANT_TOOL_PARAMS_MAP['freebayes']
        assert variant_params['tool_name'] == 'freebayes'
        find_variants_with_tool(self.alignment_group, variant_params,
                project=self.project)

        # Grab the resulting variants.
        variants = Variant.objects.filter(
                reference_genome=self.REFERENCE_GENOME)

        self._freebayes_checker(variants)

    def _freebayes_checker(self, variants):

        # There should be a Variant object for each record.
        self.assertEqual(self.EXPECTED_NUM_VARIANTS, len(variants))

        # Check that each variant is accounted for.
        self.assertEqual(set(self.EXPECTED_VARIANT_POSITIONS),
                set([v.position for v in variants]))

    def test_haploid_freebayes(self):
        """
        Perform the freebayes_test with haploid calling enabled.
        """
        self._create_alignment(haploid=True)
        self._freebayes_tester(haploid=True)

    def test_diploid_freebayes(self):
        """
        Perform the freebayes_test with haploid calling disabled.
        """
        self._create_alignment(haploid=False)
        self._freebayes_tester(haploid=False)

    def test_merged_freebayes(self):
        """
        Perform freebayes on many regions separately, and then merge.
        Compare to whole bam vcf generated in test_diploid_freebayes.
        """

        self._create_alignment(haploid=False)

        # determine the regions to send to each freebayes worker
        fb_regions = freebayes_regions(self.REFERENCE_GENOME, region_size=200)

        params = VARIANT_TOOL_PARAMS_MAP['freebayes']

        variant_param_list = []

        for region_num, fb_region in enumerate(fb_regions):
            this_region_num = region_num
            region_params = dict(params)
            region_params['tool_kwargs'] = {
                        'region':fb_region,
                        'region_num':this_region_num}
            variant_param_list.append(region_params)

        for variant_params in variant_param_list:
            find_variants_with_tool(self.alignment_group, variant_params,
                    project=self.project)

        merge_freebayes_parallel(self.alignment_group)

        # Grab the resulting variants.
        variants = Variant.objects.filter(
                reference_genome=self.REFERENCE_GENOME)

        self._freebayes_checker(variants)


class TestSVCallers(TestCase):

    def setUp(self):
        user = User.objects.create_user('test_username', password='password',
                email='test@example.com')
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')

        # Create a sample.
        self.experiment_sample = ExperimentSample.objects.create(
                project=self.project, label='sample1')

        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        self.alignment_group = alignment_group

        # Create the expected models.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=self.experiment_sample)
        bwa_dataset = Dataset.objects.create(
                label=Dataset.TYPE.BWA_ALIGN,
                type=Dataset.TYPE.BWA_ALIGN,
                status=Dataset.STATUS.READY)
        bwa_dataset.filesystem_location = clean_filesystem_location(
                TEST_DISC_SPLIT_BAM)
        bwa_dataset.save()

        sample_alignment.dataset_set.add(bwa_dataset)
        sample_alignment.save()

        self.bwa_dataset = bwa_dataset
        self.sample_alignment = sample_alignment

    def test_run_lumpy(self):
        fasta_ref = get_dataset_with_type(
            self.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        sample_alignments = [self.sample_alignment]

        vcf_output_dir = self.alignment_group.get_model_data_dir()

        vcf_output_filename = os.path.join(vcf_output_dir,'lumpy.vcf')

        alignment_type = 'BWA_ALIGN'

        bwa_disc_dataset = get_discordant_read_pairs(self.sample_alignment)
        bwa_sr_dataset = get_split_reads(self.sample_alignment)

        run_lumpy(fasta_ref, sample_alignments, vcf_output_dir,
                vcf_output_filename, alignment_type)

        dataset = Dataset.objects.create(
                type=Dataset.TYPE.VCF_LUMPY,
                label=Dataset.TYPE.VCF_LUMPY,
                filesystem_location=vcf_output_filename,
        )

        self.alignment_group.dataset_set.add(dataset)

        # Parse the resulting vcf, grab variant objects
        parse_alignment_group_vcf(self.alignment_group, Dataset.TYPE.VCF_LUMPY)

        # Grab the resulting variants.
        variants = Variant.objects.filter(reference_genome=self.reference_genome)

        # There should be a Variant object for each sv event.
        self.assertEqual(6, len(variants))

        # One event should be located very close to 25k
        va_positions = [v.position for v in variants]
        va_offset = [25000 - va_pos for va_pos in va_positions]
        self.assertTrue(any([v < 50 for v in va_offset]))

        # Clean up.
        remove_dataset_types = [
            Dataset.TYPE.LUMPY_INSERT_METRICS_MEAN_STDEV,
            Dataset.TYPE.LUMPY_INSERT_METRICS_HISTOGRAM
        ]
        for dataset_type in remove_dataset_types:
            dataset = get_dataset_with_type(self.sample_alignment, dataset_type)
            os.remove(dataset.get_absolute_location())
