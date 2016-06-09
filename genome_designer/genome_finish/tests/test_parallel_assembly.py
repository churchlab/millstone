import os

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from genome_finish.assembly import run_de_novo_assembly_pipeline
from genome_finish.millstone_de_novo_fns import create_de_novo_variants_set
from main.model_utils import get_dataset_with_type
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.testing_util import are_fastas_same
from pipeline.pipeline_runner import run_pipeline
from utils.import_util import add_dataset_to_entity
from utils.import_util import import_reference_genome_from_local_file
from utils.reference_genome_maker_util import generate_new_reference_genome

TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'

GF_TEST_DIR = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test')


class TestParallelAssembly(TestCase):

    def setUp(self):
        # Useful models.
        self.user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)
        self.project = Project.objects.create(
            owner=self.user.get_profile(), title='Test Project')

    def _get_dict_from_dir(self, data_dir):
        return {
            'ref_fasta': os.path.join(data_dir, 'ref.fa'),
            'fq_1': os.path.join(data_dir, 'reads.1.fq'),
            'fq_2': os.path.join(data_dir, 'reads.2.fq'),
            'target_fasta': os.path.join(data_dir, 'no_ins_ref.fa')
        }

    def _create_samples(self, fq_1, fq_2, num=1):

        sample_list = []
        for sample_num in range(num):
            sample = ExperimentSample.objects.create(
                    project=self.project,
                    label='test_sample_' + str(sample_num))

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

            sample_list.append(sample)

        return sample_list

    def _import_reference(self, ref_path):

        # Import reference genome
        ref_genome = import_reference_genome_from_local_file(
                self.project, 'test_ref',
                ref_path, 'fasta', move=False)
        return ref_genome

    def _align_and_assemble(self, ref_genome, sample_list):
        # Run alignment of sample to reference
        alignment_group_label = 'test_alignment'
        alignment_group, _, _ = run_pipeline(
                alignment_group_label, ref_genome, sample_list,
                perform_variant_calling=False, alignment_options={})

        # Get resulting ExperimentSampleToAlignment
        sample_align_list = ExperimentSampleToAlignment.objects.filter(
                alignment_group=alignment_group,
                experiment_sample__in=sample_list)

        # Run pipeline and wait on result
        run_de_novo_assembly_pipeline(sample_align_list)

        return alignment_group

    def _run_genome_finish_test(self, variant_set, target_fasta,
            mismatch_tolerance=0):

        self.assertTrue(variant_set.variants.exists(),
            'No placeable contigs found.')

        # Make new reference genome
        new_ref_genome_params = {'label': 'new_ref'}
        new_ref_genome = generate_new_reference_genome(
                variant_set, new_ref_genome_params)

        # Verify insertion was placed correctly
        new_ref_genome_fasta = get_dataset_with_type(
                new_ref_genome, Dataset.TYPE.REFERENCE_GENOME_FASTA
                ).get_absolute_location()

        fastas_same, indexes = are_fastas_same(
                target_fasta, new_ref_genome_fasta)

        indexes_str = str(indexes) if len(indexes) < 50 else (
                str(indexes[:50]) + '...')

        self.assertTrue(len(indexes) <= mismatch_tolerance,
                'Fastas dissimilar at indexes:' +
                indexes_str)

    def test_1_sample(self):
        data_dir = os.path.join(GF_TEST_DIR, 'small_mg1655_data/1kb_ins')
        data_dict = self._get_dict_from_dir(data_dir)

        ref_genome = self._import_reference(data_dict['ref_fasta'])
        sample_list = self._create_samples(
                data_dict['fq_1'], data_dict['fq_2'])

        alignment_group = self._align_and_assemble(ref_genome, sample_list)

        # Get set of de novo variants
        variant_set = create_de_novo_variants_set(alignment_group,
                'de_novo_variants')

        self.assertEqual(variant_set.variants.count(), 1)
        variant = variant_set.variants.all()[0]

        variant_alts = variant.variantalternate_set.all()
        self.assertEqual(len(variant_alts), 1)
        variant_alt = variant_alts[0]

        evidence_samples = [ve.experiment_sample for ve in
                variant_alt.variantevidence_set.all()]

        self.assertEqual(
                sorted(evidence_samples, key=lambda x: x.uid),
                sorted(sample_list, key=lambda x: x.uid))

    def test_2_sample(self):
        data_dir = os.path.join(GF_TEST_DIR, 'small_mg1655_data/1kb_ins')
        data_dict = self._get_dict_from_dir(data_dir)

        ref_genome = self._import_reference(data_dict['ref_fasta'])
        sample_list = self._create_samples(
                data_dict['fq_1'], data_dict['fq_2'], 2)

        alignment_group = self._align_and_assemble(ref_genome, sample_list)

        # Get set of de novo variants
        variant_set = create_de_novo_variants_set(alignment_group,
                'de_novo_variants')

        self.assertEqual(variant_set.variants.count(), 1)
        variant = variant_set.variants.all()[0]

        variant_alts = variant.variantalternate_set.all()
        self.assertEqual(len(variant_alts), 1)
        variant_alt = variant_alts[0]

        evidence_samples = [ve.experiment_sample for ve in
                variant_alt.variantevidence_set.all()]

        self.assertEqual(
                sorted(evidence_samples, key=lambda x: x.uid),
                sorted(sample_list, key=lambda x: x.uid))
