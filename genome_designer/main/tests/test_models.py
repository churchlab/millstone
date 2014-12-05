"""
Tests for models.py.
"""

import json
import os

from django.conf import settings
from django.test import TestCase

from utils.import_util import import_reference_genome_from_local_file
from utils.import_util import copy_dataset_to_entity_data_dir
from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import User
from main.models import Variant
from main.models import VariantCallerCommonData
from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.testing_util import create_common_entities
from utils import uppercase_underscore
import subprocess


TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'
TEST_PROJECT_NAME = 'testModels_project'
TEST_REF_GENOME_NAME = 'mg1655_partial'
TEST_REF_GENOME_PATH = os.path.join(settings.PWD,
    'test_data/full_vcf_test_set/mg1655_tolC_through_zupT.gb')


class TestModels(TestCase):

    def setUp(self):
        """Override.
        """
        common_entities = create_common_entities()
        self.user = common_entities['user']

    def test_delete(self):
        """Test deleting models and their associated data.

        This test was written in response to an error being thrown
        when deleting data:
        https://github.com/churchlab/genome-designer-v2/issues/219
        """
        # TODO: Add more models until we started reproducing issue #219
        # when we try to delete.
        self.user.delete()



class TestAlignmentGroup(TestCase):

    def test_get_samples(self):
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.test_project = Project.objects.create(
                title=TEST_PROJECT_NAME,
                owner=user.get_profile())
        self.test_ref_genome = ReferenceGenome.objects.create(
                project=self.test_project,
                label='blah')
        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.test_ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

        # Create a bunch of samples and relate them.
        for sample_idx in range(10):
            sample = ExperimentSample.objects.create(
                    uid=str(sample_idx),
                    project=self.test_project,
                    label='some label'
            )
            ExperimentSampleToAlignment.objects.create(
                    alignment_group=alignment_group,
                    experiment_sample=sample)

        # Test the method.
        samples = alignment_group.get_samples()
        sample_uid_set = set([sample.uid for sample in samples])
        self.assertEqual(sample_uid_set,
                set([str(x) for x in range(10)]))


class TestDataset(TestCase):

    def test_get_related_model_set(self):
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.test_project = Project.objects.create(
                title=TEST_PROJECT_NAME,
                owner=user.get_profile())
        self.test_ref_genome = ReferenceGenome.objects.create(
                project=self.test_project,
                label='blah')
        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.test_ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)
        dataset = Dataset.objects.create(
            label='the label', type=Dataset.TYPE.VCF_FREEBAYES)
        alignment_group.dataset_set.add(dataset)

        alignment_group_set = dataset.get_related_model_set()
        self.assertTrue(alignment_group in alignment_group_set.all())

    def test_dataset_compression_piping(self):
        """
        Make sure data set compression behaves correctly.
        """
        dataset = Dataset.objects.create(
                label='test_dataset',
                type=Dataset.TYPE.FASTQ1)

        GZIPPED_FASTQ_FILEPATH = os.path.join(settings.PWD, 'test_data',
                'compressed_fastq', 'sample0.simLibrary.1.fq.gz')

        dataset.filesystem_location = clean_filesystem_location(
                    GZIPPED_FASTQ_FILEPATH)

        assert dataset.is_compressed()

        process = subprocess.Popen(
                ('head '+dataset.wrap_if_compressed()+' | wc -l'),
                shell=True, executable=settings.BASH_PATH, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

        wc_output, errmsg = process.communicate()
        rc = process.returncode

        assert rc == 0, (
        "Compression process returned non-zero exit status: %s" % (
                errmsg))

        assert int(wc_output) == 10, (
                "Compression failed: %s" % (errmsg))

    def test_compress_dataset(self):
        """
        Make sure that compressing a dataset and putting a new dataset
        entry into the db works correctly.
        """
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)

        self.test_project = Project.objects.create(
            title=TEST_PROJECT_NAME,
            owner=user.get_profile())

        self.test_ref_genome = import_reference_genome_from_local_file(
            self.test_project,
            TEST_REF_GENOME_NAME,
            TEST_REF_GENOME_PATH,
            'genbank')

        dataset = get_dataset_with_type(self.test_ref_genome,
                type= Dataset.TYPE.REFERENCE_GENOME_GENBANK)

        # All the magic happens here
        compressed_dataset = dataset.make_compressed('.gz')

        # Grab the new compressed dataset through the ref genome to
        # make sure that it got added
        compressed_dataset_through_ref_genome = get_dataset_with_type(
                entity= self.test_ref_genome,
                type= Dataset.TYPE.REFERENCE_GENOME_GENBANK,
                compressed= True)
        assert compressed_dataset == compressed_dataset_through_ref_genome

    def test_dataset_strings(self):

        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)

        self.test_project = Project.objects.create(
            title=TEST_PROJECT_NAME,
            owner=user.get_profile())

        self.test_ref_genome = import_reference_genome_from_local_file(
            self.test_project,
            TEST_REF_GENOME_NAME,
            TEST_REF_GENOME_PATH,
            'genbank')

        dataset = get_dataset_with_type(self.test_ref_genome,
                type= Dataset.TYPE.REFERENCE_GENOME_GENBANK)

        self.assertEquals(
                dataset.internal_string(self.test_ref_genome),
                (str(self.test_ref_genome.uid) +
                        '_' + uppercase_underscore(Dataset.TYPE.REFERENCE_GENOME_GENBANK)))


class TestModelsStatic(TestCase):
    """Tests for static methods.
    """

    def test_clean_filesystem_location(self):
        FAKE_ABS_ROOT = '/root/of/all/evil'
        EXPECTED_CLEAN_URL = 'projects/blah'
        dirty_full_url = os.path.join(FAKE_ABS_ROOT, settings.MEDIA_ROOT,
                EXPECTED_CLEAN_URL)
        clean_location = clean_filesystem_location(dirty_full_url)
        self.assertEqual(EXPECTED_CLEAN_URL, clean_location)


class TestVariantCallerCommonData(TestCase):

    def test_json_data_field(self):
        """Tests the data field which uses the Postgresql 9.3 json type.
        """
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)

        test_project = Project.objects.create(
            title=TEST_PROJECT_NAME,
            owner=user.get_profile())

        reference_genome = ReferenceGenome.objects.create(
            project=test_project,
            label='ref1')

        chromosome = Chromosome.objects.create(
            reference_genome=reference_genome,
            label='Chromosome',
            num_bases=9001)

        variant = Variant.objects.create(
            reference_genome=reference_genome,
            type='UNKNOWN',
            chromosome=chromosome,
            position=100,
            ref_value='A'
        )

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=reference_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

        raw_data_dict = {
            'key1': 'val1',
            'key2': 'val2',
        }

        # Test storing as dictionary.
        vccd = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset_id=1,
            alignment_group=alignment_group,
            data=raw_data_dict
        )
        vccd_lookup = VariantCallerCommonData.objects.get(
            id=vccd.id)
        self.assertEquals(raw_data_dict, vccd_lookup.data)

        # Test storing as string.
        vccd = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset_id=1,
            alignment_group=alignment_group,
            data=json.dumps(raw_data_dict)
        )
        vccd_lookup = VariantCallerCommonData.objects.get(
            id=vccd.id)
        self.assertEquals(raw_data_dict, vccd_lookup.data)

        # Test blank value.
        vccd = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset_id=1,
            alignment_group=alignment_group,
        )
        self.assertEquals(0, len(vccd.data))

        # Test assigning after initial create.
        vccd = VariantCallerCommonData.objects.create(
            variant=variant,
            source_dataset_id=1,
            alignment_group=alignment_group,
        )
        vccd.data=json.dumps(raw_data_dict)
        vccd.save()
        vccd_lookup = VariantCallerCommonData.objects.get(
            id=vccd.id)
        self.assertEquals(raw_data_dict, vccd_lookup.data)


class TestExperimentSample(TestCase):

    def setUp(self):
        """Override.
        """
        self.common_entities = create_common_entities()
        self.ref_genome = self.common_entities['reference_genome']

    def test_data_dir_create_and_delete(self):
        """Make sure data directory gets deleted.
        """
        es = ExperimentSample.objects.create(
            project=self.ref_genome.project, label='test_es')

        es_data_dir = es.get_model_data_dir()

        self.assertTrue(os.path.exists(es_data_dir))

        es.delete()

        self.assertFalse(os.path.exists(es_data_dir))

    def test_add_child(self):
        """
        Make sure parent/child relationships work.
        """

        assert len(self.common_entities['sample_1'].get_children()) == 0
        assert len(self.common_entities['sample_1'].get_parents()) == 0

        self.common_entities['sample_1'].add_child(
            self.common_entities['sample_2'])

        assert len(self.common_entities['sample_1'].get_children()) == 1
        assert len(self.common_entities['sample_2'].get_parents()) == 1

        assert(self.common_entities['sample_1'].get_children()[0].uid == (
        self.common_entities['sample_2'].uid))


class TestChromosome(TestCase):

    def test_multiple_chromosome_dataset_import(self):
        user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

        project = Project.objects.create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())

        test_yeast_genome = ReferenceGenome.objects.create(
            project=project,
            label='superbrewer2000')

        test_dataset_path = os.path.join(settings.PWD, 'test_data/yeast_chrom_jkl.fasta')
        dataset_path = copy_dataset_to_entity_data_dir(test_yeast_genome, test_dataset_path)

        test_chroms_dataset  = Dataset.objects.create(
            label='jkl_chroms',
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA,
            filesystem_location=clean_filesystem_location(dataset_path))

        test_yeast_genome.dataset_set.add(test_chroms_dataset)

        # Assert correct number of chromosomes
        assert(test_yeast_genome.num_chromosomes == 3)

        # Assert correct number of bases
        assert(test_yeast_genome.num_bases == sum([chrom.num_bases for chrom in
                Chromosome.objects.filter(reference_genome=test_yeast_genome)]))

        # Assert correct chromosome labels
        expected_chrom_names = [
                'gi|448092123|ref|NC_020215.1|',
                'gi|448096713|ref|NC_020216.1|',
                'gi|448100869|ref|NC_020217.1|']

        assert([chrom.label for chrom in Chromosome.objects.filter(reference_genome=test_yeast_genome)] == expected_chrom_names)
