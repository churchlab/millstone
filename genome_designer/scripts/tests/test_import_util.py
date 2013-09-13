"""
Tests for import_util.py.
"""

import os

from django.contrib.auth.models import User
from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import VariantSet
from scripts.import_util import import_samples_from_targets_file
from scripts.import_util import import_variant_set_from_vcf
from settings import PWD as GD_ROOT_PATH


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'


class TestImportSamplesFromTargetsFile(TestCase):
    """Tests for scripts.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        # Test models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        test_project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')

    def test_import_samples(self):
        """Tests importing samples from a template file.
        """
        TARGETS_TEMPLATE_FILEPATH = os.path.join(GD_ROOT_PATH, 'main',
                'templates', 'sample_list_targets_template.tsv')

        NUM_SAMPLES_IN_TEMPLATE = 10

        # Grab any project from the database.
        project = Project.objects.all()[0]

        num_experiment_samples_before = len(ExperimentSample.objects.all())
        num_datasets_before = len(Dataset.objects.all())

        # Perform the import.
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(project,
                    UploadedFile(targets_file_fh))

        num_experiment_samples_after = len(ExperimentSample.objects.all())
        num_datasets_after = len(Dataset.objects.all())

        # Make sure the right amount of models were added.
        self.assertEqual(NUM_SAMPLES_IN_TEMPLATE,
                num_experiment_samples_after - num_experiment_samples_before)
        self.assertEqual(2 * NUM_SAMPLES_IN_TEMPLATE,
                num_datasets_after - num_datasets_before)

        # TODO: Check the filepaths as well.

class TestImportVariantSetFromVCFFile(TestCase):
    """Tests for scripts.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        # Test models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        test_project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        ref_genome = ReferenceGenome.objects.create(project=test_project,
                label='refgenome', num_chromosomes=1, num_bases=1000)


    def test_import_variant_set(self):
        """Tests importing variant sets from a pared-down vcf file
        containing only chromosome, position info, etc.
        """

        VARIANT_SET_VCF_FILEPATH = os.path.join(GD_ROOT_PATH,
                'test_data', 'fake_genome_and_reads',
                'test_genome_variant_set.vcf')

        NUM_VARIANTS_IN_SET = 20

        VARIANT_SET_NAME = 'Test Set'

        # Grab the test project / ref_genome from the database.
        # TODO: when we start checking variant integrity, make sure this is
        # the right project/reference genome!
        project = Project.objects.all()[0]
        REF_GENOME_UID = project.referencegenome_set.all()[0].uid

        ref_genome = project.referencegenome_set.all()[0]

        import_variant_set_from_vcf(ref_genome, VARIANT_SET_NAME,
                VARIANT_SET_VCF_FILEPATH)

        new_variant_set = VariantSet.objects.get(
                reference_genome=ref_genome,
                label=VARIANT_SET_NAME)

        self.assertEqual(NUM_VARIANTS_IN_SET,
                len(new_variant_set.variants.all()))

