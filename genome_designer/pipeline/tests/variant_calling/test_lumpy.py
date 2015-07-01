"""
Tests for lumpy.py.
"""

import os
import tempfile

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
from main.model_utils import clean_filesystem_location
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_split_reads
from pipeline.variant_calling import find_variants_with_tool
from pipeline.variant_calling.lumpy import filter_lumpy_vcf
from pipeline.variant_calling.lumpy import run_lumpy
from pipeline.variant_calling import TOOL_LUMPY
from pipeline.variant_calling import VARIANT_TOOL_PARAMS_MAP
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import import_reference_genome_from_local_file
from variants.vcf_parser import parse_alignment_group_vcf


TEST_DATA_DIR = os.path.join(settings.PWD, 'test_data')

TEST_FASTA = os.path.join(TEST_DATA_DIR, 'fake_genome_and_reads',
        'test_genome.fa')

TEST_DISC_SPLIT_BAM = os.path.join(settings.PWD, 'test_data',
        'discordant_split_reads', 'bwa_align.bam')

TEST_LUMPY_VCF = os.path.join(settings.PWD, 'test_data', 'pipeline',
        'variant_calling', 'lumpy.vcf')


class TestLumpy(TestCase):

    def test_run_lumpy(self):
        TEST_SAMPLE_UID = '8c57e7b9'

        user = User.objects.create_user('test_username', password='password',
                email='test@example.com')
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')

        # Create a sample.
        self.experiment_sample = ExperimentSample.objects.create(
                uid=TEST_SAMPLE_UID, project=self.project, label='sample1')

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

        fasta_ref = get_dataset_with_type(
            self.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        sample_alignments = [self.sample_alignment]

        vcf_output_dir = self.alignment_group.get_model_data_dir()

        vcf_output_filename = os.path.join(vcf_output_dir, 'lumpy.vcf')

        alignment_type = 'BWA_ALIGN'

        # NOTE: Running these functions but not checking results.
        get_discordant_read_pairs(self.sample_alignment)
        get_split_reads(self.sample_alignment)

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
        self.assertEqual(2, len(variants))

        # One event should be located very close to 25k
        va_positions = [v.position for v in variants]
        va_offset = [25000 - va_pos for va_pos in va_positions]
        self.assertTrue(any([v < 50 for v in va_offset]))


    def test_run_lumpy__deletion(self):
        """Tests running Lumpy on data that should have a deletion.
        """
        TEST_SAMPLE_UID = '38d786f2'

        user = User.objects.create_user('test_username_sv', password='password',
                email='test@example.com')

        # Grab a project.
        self.project = Project.objects.create(title='test project',
                owner=user.get_profile())

        # Use genome with deletion from our sv testing repo:
        # https://github.com/churchlab/structural-variants-testing
        DELETION_TEST_DATA_DIR = os.path.join(TEST_DATA_DIR,
                'sv_testing', 'deletion_bd5a1123')
        REF = os.path.join(DELETION_TEST_DATA_DIR, 'small_ref.fa')
        FASTQ1 = os.path.join(DELETION_TEST_DATA_DIR, 'deletion_bd5a1123.1.fq')
        FASTQ2 = os.path.join(DELETION_TEST_DATA_DIR, 'deletion_bd5a1123.2.fq')
        BWA_ALIGNMENT = os.path.join(DELETION_TEST_DATA_DIR,
                'deletion_bd5a1123.bam')

        # Create Datasets / import data.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', REF, 'fasta')
        self.experiment_sample = ExperimentSample.objects.create(
                project=self.project, label='sample1')
        copy_and_add_dataset_source(self.experiment_sample,
                Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ1, FASTQ1)
        copy_and_add_dataset_source(self.experiment_sample,
                Dataset.TYPE.FASTQ2, Dataset.TYPE.FASTQ2, FASTQ2)

        # Create an alignment that's already complete, so we can focus on
        # testing variant calling only.
        self.alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        sample_1 = ExperimentSample.objects.create(
                uid=TEST_SAMPLE_UID,
                project=self.project,
                label='sample1')

        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=self.alignment_group,
                experiment_sample=sample_1)
        copy_and_add_dataset_source(sample_alignment, Dataset.TYPE.BWA_ALIGN,
                Dataset.TYPE.BWA_ALIGN, BWA_ALIGNMENT)

        # Run lumpy.
        lumpy_params = dict(VARIANT_TOOL_PARAMS_MAP[TOOL_LUMPY])
        lumpy_params['tool_kwargs'] = {
            'region_num': sample_alignment.id,
            'sample_alignments': [sample_alignment]
        }
        find_variants_with_tool(
                self.alignment_group, lumpy_params, project=self.project)

        # Grab the resulting variants.
        variants = Variant.objects.filter(
                reference_genome=self.reference_genome)

        # Verify that we have the expected deletion around position 10000 of
        # size 1000.
        self.assertEqual(1, len(variants))
        v = variants[0]

        # start position
        self.assertTrue(9950 < v.position < 10050)

        # size
        vccd = v.variantcallercommondata_set.all()[0]
        size = vccd.data['INFO_END'] - v.position
        self.assertTrue(900 < size < 1100)

        # TODO: Check SV type.

    def test_filter_vcf(self):
        """Tests filtering out noisy values from vcf.
        """
        # We want to make sure two conditions are satisfied after filtering:
        #     1) Only records that pass rules are preserved.
        #     2) Vcf file header is preserved.

        # Read original vcf.
        with open(TEST_LUMPY_VCF) as fh:
            orig_vcf_reader = vcf.Reader(fh)
            orig_vcf_header_lines = orig_vcf_reader._header_lines
            orig_record_count = _count_records_in_vcf(fh)

        # Make sure all records present.
        self.assertTrue(len(orig_vcf_header_lines) > 1)
        self.assertEqual(4, orig_record_count)

        # Create new filtered vcf.
        NEW_VCF_FILE = tempfile.NamedTemporaryFile()
        NEW_VCF_FILE_PATH = NEW_VCF_FILE.name
        filter_lumpy_vcf(TEST_LUMPY_VCF, NEW_VCF_FILE_PATH)

        # Read new vcf.
        with open(NEW_VCF_FILE_PATH) as new_vcf_fh:
            new_vcf_reader = vcf.Reader(new_vcf_fh)
            new_vcf_header_lines = new_vcf_reader._header_lines
            new_record_count = _count_records_in_vcf(new_vcf_reader)

        # 1) Make sure header is preserved.
        self.assertEqual(set(orig_vcf_header_lines), set(new_vcf_header_lines))

        # 2) Make sure only single record.
        self.assertEqual(1, new_record_count)


def _count_records_in_vcf(vcf_reader):
    record_count = 0
    for record in vcf_reader:
        record_count += 1
    return record_count
