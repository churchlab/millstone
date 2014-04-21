"""
Tests for import_util.py.
"""

import os

from django.conf import settings
from django.contrib.auth.models import User
from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

from main.exceptions import ValidationException
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantSet
from main.testing_util import create_common_entities
from utils.import_util import DataImportError
from utils.import_util import create_sample_models_for_eventual_upload
from utils.import_util import import_reference_genome_from_local_file
from utils.import_util import import_samples_from_targets_file
from utils.import_util import import_variant_set_from_vcf
from utils.import_util import import_reference_genome_from_ncbi
from utils import internet_on

TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_DATA_ROOT = os.path.join(settings.PWD, 'test_data')
IMPORT_UTIL_TEST_DATA = os.path.join(TEST_DATA_ROOT, 'import_util_test_data')


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
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(self.project,
                    UploadedFile(targets_file_fh))

        new_samples = ExperimentSample.objects.all()
        num_experiment_samples_after = len(new_samples)
        num_datasets_after = len(Dataset.objects.all())

        # Make sure the right amount of models were added.
        self.assertEqual(NUM_SAMPLES_IN_TEMPLATE, num_experiment_samples_after)
        self.assertEqual(2 * NUM_SAMPLES_IN_TEMPLATE, num_datasets_after)

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

    def test_import_samples__no_extra_cols(self):
        """Tests importing samples from a template file that doesn't have
        extra column data filled in.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_targets_template_no_extra_col.tsv')

        # Grab any project from the database.
        project = Project.objects.all()[0]

        # Perform the import.
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(project,
                    UploadedFile(targets_file_fh))

    def test_import_samples__extra_cols(self):
        """Tests importing samples from a template file that has
        extra columns.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_targets_extra_col.tsv')

        # Grab any project from the database.
        project = Project.objects.all()[0]

        # Perform the import.
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            samples = import_samples_from_targets_file(project,
                    UploadedFile(targets_file_fh))

        # Check that the metadata was inserted successfully.
        for s in samples:
            self.assertTrue('SAMPLE_PARENT_SAMPLES' in s.data)
            self.assertTrue('SAMPLE_GROWTH_RATE' in s.data)
            self.assertTrue('SAMPLE_CYCLE' in s.data)

    def test_import_samples__bad_input(self):
        """Input data with duplicated filenames.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(IMPORT_UTIL_TEST_DATA,
                'sample_list_targets_with_duplicates.tsv')

        with self.assertRaises(AssertionError):
            with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
                import_samples_from_targets_file(self.project,
                        UploadedFile(targets_file_fh))

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
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(self.project,
                    UploadedFile(targets_file_fh))

        new_samples = ExperimentSample.objects.all()
        num_experiment_samples_after = len(new_samples)
        num_datasets_after = len(Dataset.objects.all())

        # Make sure the right amount of models were added.
        self.assertEqual(NUM_SAMPLES_IN_TEMPLATE, num_experiment_samples_after)
        self.assertEqual(2 * NUM_SAMPLES_IN_TEMPLATE, num_datasets_after)

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
                label='refgenome', num_chromosomes=1, num_bases=1000)



    def test_import_variant_set(self):
        """Tests importing variant sets from a pared-down vcf file
        containing only chromosome, position info, etc.
        """

        VARIANT_SET_VCF_FILEPATH = os.path.join(settings.PWD,
                'test_data', 'fake_genome_and_reads',
                'test_genome_variant_set.vcf')

        NUM_VARIANTS_IN_SET = 20

        VARIANT_SET_NAME = 'Test Set'

        import_variant_set_from_vcf(
                self.ref_genome, VARIANT_SET_NAME,
                VARIANT_SET_VCF_FILEPATH)

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
        self.assertEqual(set(['C','G']), set(v_553.get_alternates()))
