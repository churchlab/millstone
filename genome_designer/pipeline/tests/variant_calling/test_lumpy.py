"""
Tests for lumpy.py.
"""

import os
import tempfile

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from django.conf import settings
from django.test import TestCase
import vcf

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import get_dataset_with_type
from main.models import Variant
from main.model_utils import clean_filesystem_location
from main.testing_util import create_common_entities
from main.testing_util import create_sample_and_alignment
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_split_reads
from pipeline.variant_calling import find_variants_with_tool
from pipeline.variant_calling.lumpy import merge_lumpy_vcf
from pipeline.variant_calling.lumpy import run_lumpy
from pipeline.variant_calling.lumpy import process_vcf_post_l_merge
from pipeline.variant_calling import TOOL_LUMPY
from pipeline.variant_calling import VARIANT_TOOL_PARAMS_MAP
from utils.import_util import import_reference_genome_from_local_file
from variants.vcf_parser import parse_alignment_group_vcf
from variants.vcf_parser import SV_REF_VALUE


TEST_DATA_DIR = os.path.join(settings.PWD, 'test_data')

TEST_FASTA = os.path.join(TEST_DATA_DIR, 'fake_genome_and_reads',
        'test_genome.fa')

TEST_DISC_SPLIT_BAM = os.path.join(settings.PWD, 'test_data',
        'discordant_split_reads', 'bwa_align.bam')

# Test genomes generated using:
# https://github.com/churchlab/structural-variants-testing
DELETION_TEST_DATA_DIR = os.path.join(TEST_DATA_DIR,
        'sv_testing', 'deletion_bd5a1123')

DELETION_REF = os.path.join(DELETION_TEST_DATA_DIR, 'small_ref.fa')

# NOTE: Generated below.
DELETION_REF_GENBANK = os.path.join(DELETION_TEST_DATA_DIR, 'small_ref.gb')

# Uncomment/modify to create test data.
# def _create_annotated_ref_genome():
#     """Creates annotated ref genome.
#     """
#     with open(DELETION_REF) as fasta_fh:
#         seq_record = SeqIO.read(fasta_fh, 'fasta', alphabet=generic_dna)
#     seq_record.features = []
#     feature = SeqFeature(
#             FeatureLocation(9800, 10200, strand=1), type='CDS', id=1)
#     feature.qualifiers['gene'] = ['geneX']
#     seq_record.features.append(feature)
#     with open(DELETION_REF_GENBANK, 'w') as fh:
#         SeqIO.write(seq_record, fh, 'genbank')
# _create_annotated_ref_genome()

DELETION_FASTQ1 = os.path.join(DELETION_TEST_DATA_DIR, 'deletion_bd5a1123.1.fq')

DELETION_FASTQ2 = os.path.join(DELETION_TEST_DATA_DIR, 'deletion_bd5a1123.2.fq')

DELETION_SAMPLE_1_UID = '38d786f2'

DELETION_SAMPLE_1_BWA = os.path.join(DELETION_TEST_DATA_DIR,
        'deletion_bd5a1123_sample_uid_38d786f2.bam')

DELETION_SAMPLE_1_UID = 'ds1'

DELETION_SAMPLE_1_BWA = os.path.join(DELETION_TEST_DATA_DIR,
        'deletion_bd5a1123_ds1.bam')

DELETION_SAMPLE_2_UID = 'ds2'

DELETION_SAMPLE_2_BWA = os.path.join(DELETION_TEST_DATA_DIR,
        'deletion_bd5a1123_ds2.bam')

DELETION_SAMPLE_3_UID = 'ds3'

DELETION_SAMPLE_3_BWA = os.path.join(DELETION_TEST_DATA_DIR,
        'deletion_bd5a1123_ds3.bam')

DELETION_f8346a99_TEST_DATA_DIR = os.path.join(
        TEST_DATA_DIR, 'sv_testing', 'deletion_f8346a99')

DELETION_SAMPLE_4_UID = 'f8346a99'

DELETION_SAMPLE_4_BWA = os.path.join(DELETION_f8346a99_TEST_DATA_DIR,
        'deletion_f8346a99.bam')

INVERSION_TEST_DATA_DIR = os.path.join(
        TEST_DATA_DIR, 'sv_testing', 'inversion_5a996d78')

INVERSION_REF = os.path.join(INVERSION_TEST_DATA_DIR, 'small_ref.fa')

INVERSION_SAMPLE_UID = 'group'

INVERSION_SAMPLE_BWA = os.path.join(INVERSION_TEST_DATA_DIR,
        'inversion_5a996d78.bam')

L_MERGE_TEST_OUTPUT = os.path.join(
        TEST_DATA_DIR, 'sv_testing', 'l_merge_test_data', 'l_merge_output.vcf')


class TestLumpy(TestCase):

    def setUp(self):
        self.common_data = create_common_entities()
        self.project = self.common_data['project']

    def test_run_lumpy(self):
        TEST_SAMPLE_UID = '8c57e7b9'

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
        # Create Datasets / import data.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', DELETION_REF_GENBANK, 'genbank')

        # Create an alignment that's already complete, so we can focus on
        # testing variant calling only.
        self.alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        r = _create_sample_and_alignment(
                self.project, self.alignment_group, DELETION_SAMPLE_1_UID,
                DELETION_SAMPLE_1_BWA)
        sample_alignment = r['sample_alignment']

        # Run lumpy.
        lumpy_params = dict(VARIANT_TOOL_PARAMS_MAP[TOOL_LUMPY])
        lumpy_params['tool_kwargs'] = {
            'region_num': sample_alignment.uid,
            'sample_alignments': [sample_alignment]
        }
        find_variants_with_tool(
                self.alignment_group, lumpy_params, project=self.project)
        merge_lumpy_vcf(self.alignment_group)

        # Grab the resulting variants.
        variants = Variant.objects.filter(
                reference_genome=self.reference_genome)

        # Verify that we have the expected deletion around position 10000 of
        # size 1000.
        self.assertEqual(1, len(variants))
        v = variants[0]

        # position/ref
        self.assertTrue(9950 < v.position < 10050)
        self.assertEqual(SV_REF_VALUE, v.ref_value)

        vccd = v.variantcallercommondata_set.all()[0]

        # size
        size = vccd.data['INFO_END'] - v.position
        self.assertTrue(900 < size < 1100)

        va = v.variantalternate_set.all()[0]

        # Type
        self.assertEqual('DEL', va.data['INFO_SVTYPE'])

        # SnpEff data
        # TODO: Uncomment when Issue #648 is fixed.
        # https://github.com/churchlab/millstone/issues/648
        # self.assertEqual('geneX', va.data['INFO_EFF_GENE'])

    def test_run_lumpy__multiple_samples_of_same_exact_deletion(self):
        """Tests lumpy running on multiple samples.
        """
        # Create Datasets / import data.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', DELETION_REF, 'fasta')

        # Create an alignment that's already complete, so we can focus on
        # testing variant calling only.
        self.alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        r1 = _create_sample_and_alignment(
                self.project, self.alignment_group, DELETION_SAMPLE_1_UID,
                DELETION_SAMPLE_1_BWA)
        sa1 = r1['sample_alignment']

        r2 = _create_sample_and_alignment(
                self.project, self.alignment_group, DELETION_SAMPLE_2_UID,
                DELETION_SAMPLE_2_BWA)
        sa2 = r2['sample_alignment']

        r3 = _create_sample_and_alignment(
                self.project, self.alignment_group, DELETION_SAMPLE_3_UID,
                DELETION_SAMPLE_3_BWA)
        sa3 = r3['sample_alignment']

        r4 = _create_sample_and_alignment(
                self.project, self.alignment_group, DELETION_SAMPLE_4_UID,
                DELETION_SAMPLE_4_BWA)
        sa4 = r4['sample_alignment']

        # Common params for each run of lumpy.
        lumpy_params = dict(VARIANT_TOOL_PARAMS_MAP[TOOL_LUMPY])

        def _run_lumpy_for_sample_alignment(sa):
            """Helper function to run lumpy for sample alignment.
            """
            lumpy_params['tool_kwargs'] = {
                'region_num': sa.uid,
                'sample_alignments': [sa]
            }
            find_variants_with_tool(
                    self.alignment_group, lumpy_params, project=self.project)

        _run_lumpy_for_sample_alignment(sa1)
        _run_lumpy_for_sample_alignment(sa2)
        _run_lumpy_for_sample_alignment(sa3)
        _run_lumpy_for_sample_alignment(sa4)

        merge_lumpy_vcf(self.alignment_group)

        # Grab the resulting variants.
        variants = Variant.objects.filter(
                reference_genome=self.reference_genome)

        # Should have 2 events.
        self.assertEqual(2, len(variants))

    def test_run_lumpy__inversion(self):
        """Tests running Lumpy on data with single inversion.
        """
        # Create Datasets / import data.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', INVERSION_REF, 'fasta')

        # Create an alignment that's already complete, so we can focus on
        # testing variant calling only.
        self.alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        r = _create_sample_and_alignment(
                self.project, self.alignment_group, INVERSION_SAMPLE_UID,
                INVERSION_SAMPLE_BWA)
        sample_alignment = r['sample_alignment']

        # Run lumpy.
        lumpy_params = dict(VARIANT_TOOL_PARAMS_MAP[TOOL_LUMPY])
        lumpy_params['tool_kwargs'] = {
            'region_num': sample_alignment.uid,
            'sample_alignments': [sample_alignment]
        }
        find_variants_with_tool(
                self.alignment_group, lumpy_params, project=self.project)
        merge_lumpy_vcf(self.alignment_group)

        # Grab the resulting variants.
        variants = Variant.objects.filter(
                reference_genome=self.reference_genome)

        self.assertEqual(1, len(variants))

        v = variants[0]

        # position
        self.assertAlmostEqual(v.position, 30000, delta=2)

        # size
        vccd = v.variantcallercommondata_set.all()[0]
        size = vccd.data['INFO_END'] - v.position
        self.assertAlmostEqual(size, 1000, delta=10)

    def test_post_l_merge(self):
        """Tests post-processing code following l_sort/l_merge.py on outputs
        of lumpy applied on single samples.
        """
        _, processed_vcf_path = tempfile.mkstemp()
        process_vcf_post_l_merge(L_MERGE_TEST_OUTPUT, processed_vcf_path)

        MAP_EXPECTED_VARIANT_POS_TO_SAMPLE_LIST = {
            4998: ['f8346a99'],
            9999: ['ds1', 'ds2', 'ds3'],
        }

        with open(processed_vcf_path) as fh:
            # Assert expected number of sample cols.
            vcf_reader = vcf.Reader(fh)
            self.assertEqual(4, len(vcf_reader.samples))

            # Assert each sample has proper GT.
            for record in vcf_reader:
                samples_with_var = MAP_EXPECTED_VARIANT_POS_TO_SAMPLE_LIST[
                        record.POS]
                for sample_call in record.samples:
                    # error_msg = sample_
                    if sample_call.gt_nums == '1/1':
                        self.assertTrue(sample_call.sample in samples_with_var)
                    else:
                        self.assertFalse(sample_call.sample in samples_with_var)


###############################################################################
# Helper Functions
###############################################################################

def _count_records_in_vcf(vcf_reader):
    record_count = 0
    for record in vcf_reader:
        record_count += 1
    return record_count


def _create_sample_and_alignment(
        project, alignment_group, sample_uid, bwa_alignment):
    return create_sample_and_alignment(
            project, alignment_group, sample_uid, bwa_alignment)
