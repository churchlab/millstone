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
from pipeline.variant_calling.lumpy import filter_lumpy_vcf
from pipeline.variant_calling import run_lumpy
from utils.import_util import import_reference_genome_from_local_file
from variants.vcf_parser import parse_alignment_group_vcf


TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_DISC_SPLIT_BAM = os.path.join(settings.PWD, 'test_data',
        'discordant_split_reads', 'bwa_align.bam')

TEST_LUMPY_VCF = os.path.join(settings.PWD, 'test_data', 'pipeline',
        'variant_calling', 'lumpy.vcf')


class TestLumpy(TestCase):

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

        # Clean up.
        remove_dataset_types = [
            Dataset.TYPE.LUMPY_INSERT_METRICS_MEAN_STDEV,
            Dataset.TYPE.LUMPY_INSERT_METRICS_HISTOGRAM
        ]
        for dataset_type in remove_dataset_types:
            dataset = get_dataset_with_type(self.sample_alignment, dataset_type)
            os.remove(dataset.get_absolute_location())

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
