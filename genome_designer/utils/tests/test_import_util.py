"""
Tests for import_util.py.
"""

import os

from django.conf import settings
from django.contrib.auth.models import User
from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

from main.exceptions import ValidationException
from main.models import Chromosome
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantSet
from main.model_utils import get_dataset_with_type
from main.testing_util import create_common_entities
from utils.import_util import _get_fastqc_path
from utils.import_util import DataImportError
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import create_sample_models_for_eventual_upload
from utils.import_util import import_reference_genome_from_local_file
from utils.import_util import import_samples_from_targets_file
from utils.import_util import import_variant_set_from_vcf
from utils.import_util import import_reference_genome_from_ncbi
from utils.import_util import run_fastqc_on_sample_fastq
from utils import internet_on

TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_DATA_ROOT = os.path.join(settings.PWD, 'test_data')
IMPORT_UTIL_TEST_DATA = os.path.join(TEST_DATA_ROOT, 'import_util_test_data')

TEST_FASTQ1 = os.path.join(TEST_DATA_ROOT, 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')
TEST_FASTQ2 = os.path.join(TEST_DATA_ROOT, 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

TEST_FASTQ_GZ_1 = os.path.join(TEST_DATA_ROOT, 'fake_genome_and_reads',
        '6057f443', 'test_genome_8.snps.simLibrary.1.fq.gz')
TEST_FASTQ_GZ_2 = os.path.join(TEST_DATA_ROOT, 'fake_genome_and_reads',
        '6057f443', 'test_genome_8.snps.simLibrary.2.fq.gz')


TEST_FASTQ_GZ_B_1 = os.path.join(TEST_DATA_ROOT, 'fake_genome_and_reads',
        '70e343f5', 'test_genome_5.snps.simLibrary.1.fastq.gz')
TEST_FASTQ_GZ_B_2 = os.path.join(TEST_DATA_ROOT, 'fake_genome_and_reads',
        '70e343f5', 'test_genome_5.snps.simLibrary.2.fastq.gz')


class TestImportReferenceGenome(TestCase):
    """Tests importing a ReferenceGenome.
    """
    def setUp(self):
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                  email=TEST_EMAIL)
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')


    def test_import_reference_genome_from_local_file(self):
        """Tests importing reference genome.
        """
        TEST_GENBANK_FILE = os.path.join(settings.PWD,
                'test_data', 'import_util_test_data', 'mini_mg1655.genbank')

        import_reference_genome_from_local_file(self.project, 'a label',
                TEST_GENBANK_FILE, 'genbank')


    def test_import_reference_genome_from_local_file__fail_if_no_seq(self):
        """Should fail if no sequence in file.
        """
        TEST_GENBANK_FILE__NO_SEQ = os.path.join(settings.PWD,
                'test_data', 'import_util_test_data', 'mg1655_no_seq.genbank')

        with self.assertRaises(DataImportError):
            import_reference_genome_from_local_file(self.project, 'a label',
                    TEST_GENBANK_FILE__NO_SEQ, 'genbank')


class TestImportRefGenomeFromEntrez(TestCase):
    """Test grabbing a genome from teh interweb. Requires internet access.
    """

    def setUp(self):
        # Test models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.test_project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')

    def test_import_entrez(self):
        TEST_RECORD_ID = "6273291"
        TEST_RECORD_LABEL = "testRecord"
        if internet_on():
            import_reference_genome_from_ncbi(
                    self.test_project, TEST_RECORD_LABEL, TEST_RECORD_ID, 
                    'genbank')


class TestImportSamplesFromTargetsFile(TestCase):
    """Tests for util.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        # Test models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')

    def test_import_samples(self):
        """Tests importing samples from a template file.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(settings.PWD, 'main',
                'templates', 'sample_list_targets_template.tsv')

        NUM_SAMPLES_IN_TEMPLATE = 10

        # This test is written assuming there are no other ExperimentSamples,
        # perhaps introduced in setUp(). Check that assumption here.
        num_experiment_samples_before = len(ExperimentSample.objects.all())
        self.assertEqual(0, num_experiment_samples_before)
        num_datasets_before = len(Dataset.objects.all())
        self.assertEqual(0, num_datasets_before)

        # Perform the import.
        OPTIONS = {'skip_fastqc': True}
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(self.project,
                    UploadedFile(targets_file_fh), OPTIONS)

        new_samples = ExperimentSample.objects.all()

        # Make sure the right number of models were added.
        self.assertEqual(NUM_SAMPLES_IN_TEMPLATE, len(new_samples))
        # num_datasets_after: 2 * fastq plus
        self.assertEqual(2 * NUM_SAMPLES_IN_TEMPLATE, Dataset.objects.count())

        # Make sure all have READY status.
        observed_status_type_set = set([
                d.status for d in Dataset.objects.all()])
        self.assertEqual(1, len(observed_status_type_set))
        self.assertTrue(Dataset.STATUS.READY in observed_status_type_set)

        # Check the specifics of the Datasets, especially filepaths.
        for sample in new_samples:
            fwd_reads_dataset = sample.dataset_set.get(
                    type=Dataset.TYPE.FASTQ1)
            rev_reads_dataset = sample.dataset_set.get(
                    type=Dataset.TYPE.FASTQ2)
            self.assertNotEqual(fwd_reads_dataset.filesystem_location,
                    rev_reads_dataset.filesystem_location,
                    "Must have different filesystem locations.")

    def test_import_samples__no_extra_cols(self):
        """Tests importing samples from a template file that doesn't have
        extra column data filled in.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_targets_template_no_extra_col.tsv')

        # Grab any project from the database.
        project = Project.objects.all()[0]

        # Perform the import.
        OPTIONS = {'skip_fastqc': True}
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(project,
                    UploadedFile(targets_file_fh), OPTIONS)

    def test_import_samples__extra_cols(self):
        """Tests importing samples from a template file that has
        extra columns.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_targets_extra_col.tsv')

        # Grab any project from the database.
        project = Project.objects.all()[0]

        # Perform the import.
        OPTIONS = {'skip_fastqc': True}
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            samples = import_samples_from_targets_file(project,
                    UploadedFile(targets_file_fh), OPTIONS)

        # Check that the metadata was inserted successfully.
        for s in samples:
            self.assertTrue('SAMPLE_GROWTH_RATE' in s.data)
            self.assertTrue('SAMPLE_CYCLE' in s.data)
            self.assertTrue('SAMPLE_PARENTS' in s.data)

            # check that parents were updated
            if s.label == 'Test Sample 0':
                self.assertEqual(len(s.get_children()),4)

    def test_import_samples__bad_input(self):
        """Input data with duplicated filenames.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_targets_with_duplicates.tsv')

        OPTIONS = {'skip_fastqc': True}
        with self.assertRaises(AssertionError):
            with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
                import_samples_from_targets_file(self.project,
                        UploadedFile(targets_file_fh), OPTIONS)

    def test_import_samples__nonstandard_linebreaks(self):
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_mac_linebreaks.tsv')
        NUM_SAMPLES_IN_TEMPLATE = 10

        # This test is written assuming there are no other ExperimentSamples,
        # perhaps introduced in setUp(). Check that assumption here.
        num_experiment_samples_before = len(ExperimentSample.objects.all())
        self.assertEqual(0, num_experiment_samples_before)
        num_datasets_before = len(Dataset.objects.all())
        self.assertEqual(0, num_datasets_before)

        # Perform the import.
        # NOTE: Intentionally open as non-universal to simulate what would
        # be passed in through a browser upload.
        OPTIONS = {'skip_fastqc': True}
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(self.project,
                    UploadedFile(targets_file_fh), OPTIONS)

        new_samples = ExperimentSample.objects.all()

        # Make sure the right number of models were added.
        self.assertEqual(NUM_SAMPLES_IN_TEMPLATE, len(new_samples))
        # num_datasets_after: 2 * fastq
        self.assertEqual(2 * NUM_SAMPLES_IN_TEMPLATE, Dataset.objects.count())

        # Make sure all have READY status.
        observed_status_type_set = set([
                d.status for d in Dataset.objects.all()])
        self.assertEqual(1, len(observed_status_type_set))
        self.assertTrue(Dataset.STATUS.READY in observed_status_type_set)

        # Check the specifics of the Datasets, especially filepaths.
        for sample in new_samples:
            fwd_reads_dataset = sample.dataset_set.get(type=Dataset.TYPE.FASTQ1)
            rev_reads_dataset = sample.dataset_set.get(type=Dataset.TYPE.FASTQ2)
            self.assertNotEqual(fwd_reads_dataset.filesystem_location,
                    rev_reads_dataset.filesystem_location,
                    "Must have different filesystem locations.")


class TestCreateSampleModelsForEventualUpload(TestCase):
    """Tests the form for indicating which samples will be uploaded.
    """
    def setUp(self):
        """Override.
        """
        common_entities = create_common_entities()
        self.project = common_entities['project']

    def test_paired_reads_upload(self):
        """Basic test.
        """
        FILLED_OUT_FORM = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_browser_upload_test_data.tsv')

        num_samples_before = ExperimentSample.objects.count()

        with open(FILLED_OUT_FORM) as targets_file:
            create_sample_models_for_eventual_upload(self.project, targets_file)

        EXPECTED_NUM_SAMPLES_ADDED = 2
        samples = ExperimentSample.objects.all()
        self.assertEqual(EXPECTED_NUM_SAMPLES_ADDED,
                len(samples) - num_samples_before)

        # Make sure the Datasets have the right status and location.
        for sample_prefix in ['sample1', 'sample2']:
            sample_label = sample_prefix + '_to_upload'
            es = ExperimentSample.objects.get(label=sample_label)
            fq1_ds = es.dataset_set.get(type=Dataset.TYPE.FASTQ1)
            fq2_ds = es.dataset_set.get(type=Dataset.TYPE.FASTQ2)
            self.assertEqual(Dataset.STATUS.AWAITING_UPLOAD, fq1_ds.status)
            self.assertEqual(Dataset.STATUS.AWAITING_UPLOAD, fq2_ds.status)
            EXPECTED_FQ1_PATH = os.path.join(settings.MEDIA_ROOT,
                    es.get_model_data_dir(), sample_prefix + '_r1.fastq')
            self.assertEqual(EXPECTED_FQ1_PATH,  fq1_ds.get_absolute_location())
            EXPECTED_FQ2_PATH = os.path.join(settings.MEDIA_ROOT,
                    es.get_model_data_dir(), sample_prefix + '_r2.fastq')
            self.assertEqual(EXPECTED_FQ2_PATH, fq2_ds.get_absolute_location())

    def test_unpaired_reads_upload(self):
        FILLED_OUT_FORM = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_browser_upload_test_data__unpaired.tsv')

        num_samples_before = ExperimentSample.objects.count()

        with open(FILLED_OUT_FORM) as targets_file:
            create_sample_models_for_eventual_upload(self.project, targets_file)

        EXPECTED_NUM_SAMPLES_ADDED = 2
        samples = ExperimentSample.objects.all()
        self.assertEqual(EXPECTED_NUM_SAMPLES_ADDED,
                len(samples) - num_samples_before)

        # Make sure the Datasets have the right status and location.
        for sample_prefix in ['sample1', 'sample2']:
            sample_label = sample_prefix + '_to_upload'
            es = ExperimentSample.objects.get(label=sample_label)

            # Check fastq1.
            fq1_ds = es.dataset_set.get(type=Dataset.TYPE.FASTQ1)
            self.assertEqual(Dataset.STATUS.AWAITING_UPLOAD, fq1_ds.status)
            EXPECTED_FQ1_PATH = os.path.join(settings.MEDIA_ROOT,
                    es.get_model_data_dir(), sample_prefix + '_r1.fastq')
            self.assertEqual(EXPECTED_FQ1_PATH,  fq1_ds.get_absolute_location())

            # No Dataset for fastq2.
            self.assertEqual(0,
                    es.dataset_set.filter(type=Dataset.TYPE.FASTQ2).count())

    def test_same_filename_for_read1_and_read2(self):
        """Tests a common error I can imagine where the user forgets to fix
        filenames.
        """
        FILLED_OUT_FORM = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_browser_upload_test_data__repeated_filename.tsv')

        with self.assertRaises(ValidationException):
            with open(FILLED_OUT_FORM) as targets_file:
                create_sample_models_for_eventual_upload(self.project,
                        targets_file)

    def test_attempt_add_existing_filename(self):
        """Tests when upload form has a row that would add a sample with
        the same filename.
        """
        FILLED_OUT_FORM = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_browser_upload_test_data.tsv')

        # First upload succeeds.
        with open(FILLED_OUT_FORM) as targets_file:
            create_sample_models_for_eventual_upload(self.project, targets_file)

        # Second upload fails due to repeated filenames.
        with self.assertRaises(ValidationException):
            with open(FILLED_OUT_FORM) as targets_file:
                create_sample_models_for_eventual_upload(self.project,
                        targets_file)

    def test_form_with_repeated_names(self):
        """Tests when upload form has a row that would add a sample with
        the same filename.
        """
        FILLED_OUT_FORM = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_browser_upload_test_data__repeated_sample_names.tsv')

        # Second upload fails due to repeated filenames.
        with self.assertRaises(ValidationException):
            with open(FILLED_OUT_FORM) as targets_file:
                create_sample_models_for_eventual_upload(self.project,
                        targets_file)


class TestImportVariantSetFromVCFFile(TestCase):
    """Tests for util.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        # Test models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        test_project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=test_project,
                label='refgenome')
        Chromosome.objects.create(
            reference_genome=self.ref_genome,
            label='Chromosome',
            num_bases=9001)

    def _assert_variants(self, vcf_filepath):
        """Common assert code used by multiple tests.
        """
        NUM_VARIANTS_IN_SET = 20

        VARIANT_SET_NAME = 'Test Set'

        import_variant_set_from_vcf(
                self.ref_genome, VARIANT_SET_NAME, vcf_filepath)

        new_variant_set = VariantSet.objects.get(
                reference_genome=self.ref_genome,
                label=VARIANT_SET_NAME)

        self.assertEqual(NUM_VARIANTS_IN_SET,
                len(new_variant_set.variants.all()))

        # Spot-check a few variants.
        v_1128 = Variant.objects.get(reference_genome=self.ref_genome,
                position=1128)
        self.assertEqual(['C'], v_1128.get_alternates())

        v_553 = Variant.objects.get(reference_genome=self.ref_genome,
                position=553)
        self.assertEqual(set(['C', 'G']), set(v_553.get_alternates()))

    def test_import_variant_set(self):
        """Tests importing variant sets from a pared-down vcf file
        containing only chromosome, position info, etc.
        """
        VARIANT_SET_VCF_FILEPATH = os.path.join(settings.PWD,
                'test_data', 'fake_genome_and_reads',
                'test_genome_variant_set.vcf')
        self._assert_variants(VARIANT_SET_VCF_FILEPATH)

    def test_import_variant_set__nonstandard_linebreaks(self):
        """Tests importing variant sets from file with nonstandard linebreaks
        as might happen with a tab-separated file saved in Excel on Mac OS X.
        """
        VARIANT_SET_VCF_FILEPATH = os.path.join(settings.PWD,
                'test_data', 'fake_genome_and_reads',
                'test_genome_variant_set__mac_linebreaks.vcf')
        self._assert_variants(VARIANT_SET_VCF_FILEPATH)


class TestFastQC(TestCase):
    """Tests running fastqc util.
    """

    def setUp(self):
        self.common_entities = create_common_entities()

    def _fastqc_test_runner(self, fastq1_location, fastq2_location):
        """Helper that takes different fastqs as source.

        This function is a test itself.
        """
        # Run FastQC
        gz_backed_sample = self.common_entities['sample_1']
        gz_fastq1_dataset = copy_and_add_dataset_source(
                gz_backed_sample, Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ1,
                fastq1_location)
        gz_fastq2_dataset = copy_and_add_dataset_source(
                gz_backed_sample, Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ2,
                fastq2_location)
        run_fastqc_on_sample_fastq(gz_backed_sample, gz_fastq1_dataset)
        run_fastqc_on_sample_fastq(gz_backed_sample, gz_fastq2_dataset,
                rev=True)

        # We expect 2 Dataset per Fastq so 4 total.
        self.assertEqual(4, Dataset.objects.count())

        # Check link matches file extension.
        FASTQC_DATASET_TYPES = [
                Dataset.TYPE.FASTQC1_HTML, Dataset.TYPE.FASTQC2_HTML]
        for fastqc_dataset_type in FASTQC_DATASET_TYPES:
            fastqc_1_dataset = get_dataset_with_type(
                    gz_backed_sample, fastqc_dataset_type)
            assert os.path.exists(fastqc_1_dataset.get_absolute_location())

    def test_fastqc(self):
        self._fastqc_test_runner(TEST_FASTQ1, TEST_FASTQ2)

    def test_fastqc_gzipped(self):
        self._fastqc_test_runner(TEST_FASTQ_GZ_1, TEST_FASTQ_GZ_2)

    def test_get_fastqc_filename(self):
        FASTQ_FILENAME = 'R2_001.fq.gz'

        fastq_dataset = Dataset.objects.create(
                label='irrelevant',
                type='also irrelevant',
                filesystem_location=FASTQ_FILENAME)

        actual_fastqc_filename = os.path.split(
                _get_fastqc_path(fastq_dataset))[1]

        EXPECTED_FASTQC_FILENAME = 'R2_001.fq_fastqc.html'
        self.assertEqual(EXPECTED_FASTQC_FILENAME, actual_fastqc_filename)

    def test_get_fastqc_filename__extension_is_fastqc(self):
        FASTQ_FILENAME = 'R2_001.fastq.gz'

        fastq_dataset = Dataset.objects.create(
                label='irrelevant',
                type='also irrelevant',
                filesystem_location=FASTQ_FILENAME)

        actual_fastqc_filename = os.path.split(
                _get_fastqc_path(fastq_dataset))[1]

        EXPECTED_FASTQC_FILENAME = 'R2_001_fastqc.html'
        self.assertEqual(EXPECTED_FASTQC_FILENAME, actual_fastqc_filename)
