"""
Tests for genome finishing features
"""
import os

from Bio import SeqIO
from django.conf import settings
from django.contrib.auth import authenticate
from django.contrib.auth.models import User
from django.http.request import HttpRequest
from django.test import Client
from django.test import TestCase

from main.model_utils import get_dataset_with_type
from main.models import Contig
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import VariantSet
import main.xhr_handlers as xhr_handlers
from pipeline.pipeline_runner import run_pipeline
from utils import are_fastas_same
from utils.import_util import add_dataset_to_entity
from utils.import_util import import_reference_genome_from_local_file
from utils.reference_genome_maker_util import generate_new_reference_genome

TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'

GF_TEST_DIR = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test')


class TestGenomeFinishMG1655(TestCase):

    def setUp(self):
        # Useful models.
        self.user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)
        self.project = Project.objects.create(
            owner=self.user.get_profile(), title='Test Project')

        # Fake web browser client used to make requests.
        self.client = Client()
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)

    def _perform_assembly(self, data_dir):

        ref_fasta = os.path.join(data_dir, 'ref.fa')
        fq_1 = os.path.join(data_dir, 'reads.1.fq')
        fq_2 = os.path.join(data_dir, 'reads.2.fq')

        # Import reference genome
        ref_genome = import_reference_genome_from_local_file(
                self.project, 'test_ref',
                ref_fasta, 'fasta', move=False)

        # Create sample model
        sample = ExperimentSample.objects.create(
                project=self.project,
                label='test_sample')

        # Add fastq datasets to sample
        add_dataset_to_entity(
                sample,
                Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1,
                filesystem_location=fq_1)
        add_dataset_to_entity(
                sample,
                Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2,
                filesystem_location=fq_2)

        # Run alignment of sample to reference
        alignment_group_label = 'test_alignment'
        sample_list = [sample]
        alignment_group, _, _ = run_pipeline(
                alignment_group_label, ref_genome, sample_list,
                perform_variant_calling=False, alignment_options={})

        # Get resulting ExperimentSampleToAlignment
        sample_align = ExperimentSampleToAlignment.objects.get(
                alignment_group=alignment_group,
                experiment_sample=sample)

        # Create HttpRequest
        request = HttpRequest()
        request_data = {
            'sampleAlignmentUid': sample_align.uid
        }
        request.GET = request_data
        request.method = 'GET'
        request.user = self.user

        # Send request
        authenticate(username=TEST_USERNAME, password=TEST_PASSWORD)
        self.assertTrue(request.user.is_authenticated())
        xhr_handlers.generate_contigs(request)

        # Retrieve contigs
        contigs = Contig.objects.filter(
                parent_reference_genome=ref_genome,
                experiment_sample_to_alignment=sample_align)

        return contigs


    def test_1kb_insertion(self):

        data_dir_1kb = os.path.join(GF_TEST_DIR, 'small_mg1655_data/1kb_ins')
        contigs = self._perform_assembly(data_dir_1kb)

        # Assert one contig was generated
        self.assertEqual(contigs.count(), 1)
        self.assertTrue(contigs[0].num_bases > 0)

        # Get VariantSet
        insertion_variant_set = VariantSet.objects.get(
                reference_genome=contigs[0].parent_reference_genome,
                label='de_novo_assembled_insertions')

        # Make new reference genome
        new_ref_genome_params = {'label': 'new_ref'}
        new_ref_genome = generate_new_reference_genome(
                insertion_variant_set, new_ref_genome_params)

        # Verify insertion was placed correctly
        target_fasta = os.path.join(data_dir_1kb, 'no_ins_ref.fa')
        new_ref_genome_fasta = get_dataset_with_type(
                new_ref_genome, Dataset.TYPE.REFERENCE_GENOME_FASTA
                ).get_absolute_location()

        self.assertTrue(are_fastas_same(target_fasta, new_ref_genome_fasta))
