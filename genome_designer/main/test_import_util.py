"""
Tests for import_util.py.
"""

import os

from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

from genome_designer.main.models import Dataset
from genome_designer.main.models import ExperimentSample
from genome_designer.main.models import Project
from genome_designer.scripts.import_util import import_samples_from_targets_file
from genome_designer.scripts.bootstrap_data import bootstrap_fake_data
from genome_designer.settings import PWD as GD_ROOT_PATH

class TestImportSamplesFromTargetsFile(TestCase):
    """Tests for scripts.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        bootstrap_fake_data()

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
