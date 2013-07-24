"""
Tests for import_util.py.
"""

import os

from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

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

        # Grab any project from the database.
        project = Project.objects.all()[0]

        # Grab the template file.
        with open(TARGETS_TEMPLATE_FILEPATH) as targets_file_fh:
            import_samples_from_targets_file(project,
                    UploadedFile(targets_file_fh))
