"""
Tests for snpeff_util.py
"""

"""
Tests for alignment_pipeline.py
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
from pipeline.variant_effects import build_snpeff
from pipeline.variant_effects import run_snpeff
from pipeline.variant_effects import get_snpeff_config_path
from pipeline.variant_effects import populate_record_eff
from pipeline.variant_effects import SNPEFF_ALT_RE
from pipeline.variant_calling.common import add_vcf_dataset
from pipeline.variant_calling.constants import TOOL_FREEBAYES
from utils.import_util import import_reference_genome_from_local_file


TEST_DIR = os.path.join(settings.PWD, 'test_data', 'genbank_aligned')

TEST_GENBANK = os.path.join(TEST_DIR, 'mg1655_tolC_through_zupT.gb')

TEST_UNANNOTATED_VCF = os.path.join(TEST_DIR, 'bwa_align_unannotated.vcf')

NO_VARIANTS_VCF = os.path.join(TEST_DIR, 'no_variants.vcf')


class TestSnpeff(TestCase):

    def setUp(self):
        user = User.objects.create_user('test_username', password='password',
                email='test@example.com')

        # Grab a project.
        self.project = Project.objects.create(title='snpeff test project',
                owner=user.get_profile())

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'snpeff test ref genome', TEST_GENBANK, 'genbank')

        # Create a new alignment group.
        self.alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Create a sample.
        self.sample_1 = ExperimentSample.objects.create(
                project=self.project,
                label='test sample 1')

        # Create relationship between alignment and sample.
        self.sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=self.alignment_group,
                experiment_sample=self.sample_1)

        # Add unannotated SNP data.
        self.vcf_dataset = Dataset.objects.create(
                type=Dataset.TYPE.VCF_FREEBAYES,
                label=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location=TEST_UNANNOTATED_VCF)
        self.alignment_group.dataset_set.add(self.vcf_dataset)

    def test_build_snpeff(self):
        """ Run the config pipeline and check that all the required files,
            the reference genome config file and the snpeff database file,
            are present.
        """
        build_snpeff(self.reference_genome)

        self.assertTrue(
                os.path.exists(self.reference_genome.get_snpeff_config_path()),
                'SnpEff config file was not found.')

        snpeff_database_file = os.path.join(
                self.reference_genome.get_snpeff_genbank_parent_dir(),
                'snpEffectPredictor.bin')
        self.assertTrue(
                os.path.exists(snpeff_database_file),
                'SnpEff annotation database was not found.')

    def test_run_snpeff(self):
        """Test running the pipeline that annotates SNPS.

        This test doesn't check the accuracy of the SNP-annotation. The test is
        intended just to run the pipeline and make sure there are no errors.
        """

        # Try running snpeff
        snpeff_vcf_filename = run_snpeff(
                self.alignment_group, TOOL_FREEBAYES)

        vcf_dataset = add_vcf_dataset(
                self.alignment_group, Dataset.TYPE.VCF_FREEBAYES_SNPEFF,
                snpeff_vcf_filename)

        # Check that the alignment group has a freebayes vcf dataset associated
        # with it.
        vcf_dataset = get_dataset_with_type(
                self.alignment_group, Dataset.TYPE.VCF_FREEBAYES_SNPEFF)
        self.assertIsNotNone(vcf_dataset,
            'SnpEff annotated vcf dataset was not found after running snpeff.')

        # Make sure the .vcf file actually exists.
        self.assertTrue(os.path.exists(vcf_dataset.get_absolute_location()),
            'SnpEff annotated vcf file was not found after running snpeff.')

        # Make sure the vcf is valid by reading it using pyvcf.
        with open(vcf_dataset.get_absolute_location()) as vcf_fh:
            reader = vcf.Reader(vcf_fh)
            record = reader.next()
            assert 'EFF' in record.INFO, (
                    'No EFF INFO field found in snpeff VCF file.')
            assert not 'ERROR' in record.INFO['EFF']
            for record in reader:
                assert not 'ERROR' in record.INFO['EFF']

    def test_populate_record_eff(self):
        """ Test the regex on a few snpeff field examples.
        """

        class FakeVCFRecord(object):
            def __init__(self):
                self.INFO = {}
                self.CHROM = 'U00096.2'
                self.POS = 100

        test_record = FakeVCFRecord()

        # single eff
        test_record.INFO['EFF'] = [''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC',
            '||CODING|b3038|1|1)'))]
        updated_test_record = populate_record_eff(test_record)

        # eff with error field
        test_record.INFO['EFF'] = [''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC',
            '||CODING|b3038|1|1|WARN_TEST|ERROR_TEST)'))]
        updated_test_record = populate_record_eff(test_record)

        # missing eff field
        test_record.INFO.pop('EFF', None)
        updated_test_record = populate_record_eff(test_record)

        self.assertEqual(updated_test_record.INFO['EFF_EFFECT'], ['ERROR'])
        self.assertTrue('SNPEFF_ERROR' in updated_test_record.INFO['EFF_ERR'][0])

        # multi-eff
        test_record.INFO['EFF'] = [''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC',
            '||CODING|b3038|1|1|WARN_TEST|ERROR_TEST),')),''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aGg|M239T|386|ygiC',
            '||CODING|b3038|1|1|ERROR_TEST|WARN_TEST)'))]
        updated_test_record = populate_record_eff(test_record)

        self.assertEqual(updated_test_record.INFO['EFF_CONTEXT'],
                ['aTg/aCg','aTg/aGg'])

        self.assertEqual(updated_test_record.INFO['EFF_GENE'],
                ['ygiC','ygiC'])

    def test_variant_effects__no_variants(self):
        """Tests that variant_effects runs without problems when input
        is empty vcf.
        """
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)
        vcf_dataset = Dataset.objects.create(
                type=Dataset.TYPE.VCF_FREEBAYES,
                label=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location=NO_VARIANTS_VCF)
        alignment_group.dataset_set.add(vcf_dataset)

        run_snpeff(alignment_group, TOOL_FREEBAYES)


class TestSnpeffNoSetup(TestCase):
    """Small tests that don't require setup.
    """

    def test_snpeff_alt_re__no_EFFECT(self):
        match = SNPEFF_ALT_RE.match(
                '(MODIFIER||||||||||G|ERROR_OUT_OF_CHROMOSOME_RANGE)')
        self.assertTrue(match)
