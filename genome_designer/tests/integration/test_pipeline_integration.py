"""Alignment pipeline integration tests.
"""

import os
import time

from django.conf import settings
from django.test import TransactionTestCase

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.testing_util import create_common_entities
from pipeline.pipeline_runner import run_pipeline
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from scripts.import_util import import_reference_genome_from_ncbi
from scripts.util import internet_on


TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')


class TestAlignmentPipeline(TransactionTestCase):

    def setUp(self):
        common_entities = create_common_entities()
        self.project = common_entities['project']
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')

        self.experiment_sample = ExperimentSample.objects.create(
                project=self.project, label='sample1')
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)

    def test_run_pipeline(self):
        """Tests running the full pipeline.
        """
        sample_list = [self.experiment_sample]

        alignment_group_obj, async_result = run_pipeline('name_placeholder',
                self.reference_genome, sample_list)

        # Block until pipeline finishes.
        while not async_result.ready():
            time.sleep(1)
        if async_result.status == 'FAILURE':
            self.fail('Async task failed.')

        # Refresh the object.
        alignment_group_obj = AlignmentGroup.objects.get(
                id=alignment_group_obj.id)

        # Verify the AlignmentGroup object is created.
        self.assertEqual(1,
                len(alignment_group_obj.experimentsampletoalignment_set.all()))
        self.assertEqual(AlignmentGroup.STATUS.COMPLETED,
                alignment_group_obj.status)

        # Make sure the initial JBrowse config has been created.
        jbrowse_dir = self.reference_genome.get_jbrowse_directory_path()
        self.assertTrue(os.path.exists(jbrowse_dir))
        self.assertTrue(os.path.exists(os.path.join(jbrowse_dir,
                'indiv_tracks')))
