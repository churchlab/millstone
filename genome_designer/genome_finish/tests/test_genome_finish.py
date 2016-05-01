"""
Tests for genome finishing features
"""
import os
import re

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from genome_finish.assembly import evaluate_contigs
from genome_finish.assembly import parse_variants_from_vcf
from genome_finish.assembly import run_de_novo_assembly_pipeline
from genome_finish.millstone_de_novo_fns import create_de_novo_variants_set
from main.model_utils import get_dataset_with_type
from main.models import AlignmentGroup
from main.models import Contig
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.models import VariantSet
from main.testing_util import are_fastas_same
from pipeline.pipeline_runner import run_pipeline
from utils.import_util import add_dataset_to_entity
from utils.import_util import import_reference_genome_from_local_file
from utils.reference_genome_maker_util import generate_new_reference_genome
from variants.variant_sets import update_variant_in_set_memberships

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

    def populate_dict_from_dir(self, data_dir):
        return {
            'ref_fasta': os.path.join(data_dir, 'ref.fa'),
            'fq_1': os.path.join(data_dir, 'reads.1.fq'),
            'fq_2': os.path.join(data_dir, 'reads.2.fq'),
            'target_fasta': os.path.join(data_dir, 'no_ins_ref.fa')
        }


    def _perform_assembly(self, data_dict):

        ref_fasta = data_dict['ref_fasta']
        fq_1 = data_dict['fq_1']
        fq_2 = data_dict['fq_2']

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

        # Run pipeline and wait on result
        async_result = run_de_novo_assembly_pipeline([sample_align])
        async_result.get()

        # Retrieve contigs
        contigs = Contig.objects.filter(
                parent_reference_genome=ref_genome,
                experiment_sample_to_alignment=sample_align)

        return contigs

    def _run_genome_finish_test(self, data_dict, mismatch_tolerance=0):

        contigs = self._perform_assembly(data_dict)

        # Assert contigs were generated
        self.assertTrue(contigs.count() > 0)
        self.assertTrue(contigs[0].num_bases > 0)

        ag = contigs[0].experiment_sample_to_alignment.alignment_group

        # Get set of de novo variants
        variant_set = create_de_novo_variants_set(ag, 'de_novo_variants')

        contigs_found_error_str = (str(len(contigs)) + ' found with lengths:' +
            ', '.join([str(c.num_bases) for c in contigs]))

        self.assertTrue(variant_set.variants.exists(),
            'No placeable contigs found.  ' +
            contigs_found_error_str)

        # Make new reference genome
        new_ref_genome_params = {'label': 'new_ref'}
        new_ref_genome = generate_new_reference_genome(
                variant_set, new_ref_genome_params)

        # Verify insertion was placed correctly
        target_fasta = data_dict['target_fasta']
        new_ref_genome_fasta = get_dataset_with_type(
                new_ref_genome, Dataset.TYPE.REFERENCE_GENOME_FASTA
                ).get_absolute_location()

        fastas_same, indexes = are_fastas_same(
                target_fasta, new_ref_genome_fasta)

        indexes_str = str(indexes) if len(indexes) < 50 else (
                str(indexes[:50]) + '...')

        self.assertTrue(len(indexes) <= mismatch_tolerance,
                'Fastas dissimilar at indexes:' +
                indexes_str + '\n' +
                contigs_found_error_str)

    def test_1kb_insertion(self):
        data_dir = os.path.join(GF_TEST_DIR, 'small_mg1655_data/1kb_ins')
        self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    def test_1kb_insertion_cov_80(self):
        data_dir = os.path.join(GF_TEST_DIR,
                'small_mg1655_data/1kb_ins_cov_80')
        self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    # # Due to adapter overlap issues this test currently fails
    # def test_1kb_insertion_cov_40(self):
    #     data_dir = os.path.join(GF_TEST_DIR,
    #             'small_mg1655_data/1kb_ins_cov_40')
    #     self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    def test_1kb_insertion_del_30(self):
        data_dir = os.path.join(GF_TEST_DIR,
                'small_mg1655_data/1kb_ins_del_30')
        self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    def test_1kb_insertion_del_1000(self):
        data_dir = os.path.join(GF_TEST_DIR,
                'small_mg1655_data/1kb_ins_del_1000')
        self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    def test_gzip_reads(self):
        data_dir = os.path.join(GF_TEST_DIR, 'gz_reads')
        data_dict = {k: os.path.join(data_dir, f) for k, f in
                [('ref_fasta', 'del_genome.fa'), ('target_fasta', 'genome.fa'),
                ('fq_1', 'simngs_reads.1.fq.gz'),
                ('fq_2', 'simngs_reads.2.fq.gz')]
        }

        self._run_genome_finish_test(data_dict)

    # # Due to issues with false positive translocations, this test
    # # curently fails
    # def test_translocation(self):
    #     data_dir = os.path.join(GF_TEST_DIR, 'random_seq_data',
    #         'is_element_reads')
    #     self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    def test_translocation_far_1(self):
        data_dir = os.path.join(GF_TEST_DIR, 'random_seq_data',
                'is_element_far_1')
        self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))

    def test_translocation_far_long(self):
        data_dir = os.path.join(GF_TEST_DIR, 'random_seq_data',
                'is_element_far_long')
        self._run_genome_finish_test(self.populate_dict_from_dir(data_dir))


class TestGraphWalk(TestCase):

    def setUp(self):
        # Useful models.
        self.user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)
        self.project = Project.objects.create(
            owner=self.user.get_profile(), title='Test Project')

    def _make_dummy_models(self):

        self.reference_genome = ReferenceGenome.objects.create(
                project=self.project,
                label='test_reference_genome')


        alignment_group = AlignmentGroup.objects.create(
                reference_genome=self.reference_genome)

        sample = ExperimentSample.objects.create(
                project=self.project,
                label='test_sample')

        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=sample)

        return {'sample_alignment': sample_alignment,
                'reference_genome': self.reference_genome,
                'alignment_group': alignment_group}

    def _run_contig_walk_test(self, test_dir):

        ref_fasta = os.path.join(test_dir, 'ref.fa')
        self.target_fasta = os.path.join(test_dir, 'target.fa')

        contig_fasta_list = filter(
                lambda x: re.match(r'contig_\d+\.fa', x),
                os.listdir(test_dir))
        contig_fasta_list = [os.path.join(test_dir, filename) for
                filename in contig_fasta_list]

        dummy_models = self._make_dummy_models()
        reference_genome = dummy_models['reference_genome']
        sample_alignment = dummy_models['sample_alignment']
        alignment_group = dummy_models['alignment_group']

        add_dataset_to_entity(
                    reference_genome,
                    Dataset.TYPE.REFERENCE_GENOME_FASTA,
                    Dataset.TYPE.REFERENCE_GENOME_FASTA,
                    filesystem_location=ref_fasta)


        ref_genbank = os.path.join(test_dir, 'ref.gb')
        if os.path.exists(ref_genbank):
            add_dataset_to_entity(
                reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_GENBANK,
                Dataset.TYPE.REFERENCE_GENOME_GENBANK,
                ref_genbank)

        # Make data_dir directory to house genome_finishing files
        assembly_dir = os.path.join(
                sample_alignment.get_model_data_dir(),
                'assembly')

        # Make assembly directory
        os.mkdir(assembly_dir)

        data_dir = os.path.join(assembly_dir, '0')
        os.mkdir(data_dir)

        # Create contigs
        contig_list = []
        for i, contig_fasta in enumerate(contig_fasta_list):
            contig = Contig.objects.create(
                parent_reference_genome=reference_genome,
                experiment_sample_to_alignment=sample_alignment,
                label='test_contig_' + str(i))
            add_dataset_to_entity(
                    contig,
                    Dataset.TYPE.REFERENCE_GENOME_FASTA,
                    Dataset.TYPE.REFERENCE_GENOME_FASTA,
                    filesystem_location=contig_fasta)
            contig.metadata['assembly_dir'] = data_dir
            contig.metadata['node_number'] = i
            contig_list.append(contig)

        # Place contigs and create variants
        evaluate_contigs(contig_list,
                skip_extracted_read_alignment=True,
                use_read_alignment=False)
        parse_variants_from_vcf(sample_alignment)

        # Get set of de novo variants
        variant_set = create_de_novo_variants_set(
                alignment_group, 'de_novo_variants')

        for v in variant_set.variants.all():
            alts = v.get_alternates()
            assert len(alts) == 1
            alt = alts[0]
            print '\npos:%s\nref: %dbp :%s\nalt: %dbp :%s\n' % (
                    v.position,
                    len(v.ref_value), v.ref_value,
                    len(alt), alt)

        return variant_set

    def _assert_variants_make_target(self, variant_set):
        self.assertTrue(variant_set.variants.exists())
        # self.assertEqual(len(variant_set.variants.all()), 1)

        if len(variant_set.variants.all()) == 1:
            # Make new reference genome
            new_ref_genome_params = {'label': 'new_ref'}
            new_ref_genome = generate_new_reference_genome(
                    variant_set, new_ref_genome_params)

            # Verify insertion was placed correctly
            new_ref_genome_fasta = get_dataset_with_type(
                    new_ref_genome, Dataset.TYPE.REFERENCE_GENOME_FASTA
                    ).get_absolute_location()

            fastas_same, indexes = are_fastas_same(
                    self.target_fasta, new_ref_genome_fasta)

            self.assertTrue(fastas_same)

        else:
            print 'Multiple variants'
            incorrect_variants = []
            for i, v in enumerate(variant_set.variants.all()):

                # Make Variant set for single variant
                variant_set = VariantSet.objects.create(
                        reference_genome=self.reference_genome,
                        label='multiple_variant_%d' % i)

                update_variant_in_set_memberships(
                    self.reference_genome,
                    [v.uid],
                    'add',
                    variant_set.uid)

                new_ref_genome_params = {'label': 'new_ref'}

                new_ref_genome = generate_new_reference_genome(
                        variant_set, new_ref_genome_params)

                # Verify insertion was placed correctly
                new_ref_genome_fasta = get_dataset_with_type(
                        new_ref_genome, Dataset.TYPE.REFERENCE_GENOME_FASTA
                        ).get_absolute_location()

                fastas_same, indexes = are_fastas_same(
                        self.target_fasta, new_ref_genome_fasta)

                if not fastas_same:
                    incorrect_variants.append((i, v))

            if incorrect_variants:
                print 'Variant %d/%d was incorrect' * len(
                        incorrect_variants) * tuple(zip(
                                [t[0] for t in incorrect_variants],
                                len(incorrect_variants) *
                                len(variant_set.variants.all())))
            self.assertTrue(not incorrect_variants, '%d incorrect variants' % (
                    len(incorrect_variants) /
                    len(variant_set.variants.all())))

    def test_deletion(self):
        test_dir = os.path.join(GF_TEST_DIR, 'random_seq_data',
                'deletion')
        variant_set = self._run_contig_walk_test(test_dir)
        self._assert_variants_make_target(variant_set)

    def test_homology_flanked_deletion(self):
        test_dir = os.path.join(GF_TEST_DIR, 'random_seq_data',
                'homology_flanked_deletion')
        variant_set = self._run_contig_walk_test(test_dir)
        self._assert_variants_make_target(variant_set)

    # def test_annotated_mobile_element(self):
    #     """Test for discovery of mobile element insertion in the
    #     case where there is an integrated copy of the mobile element in the
    #     genome annotated in the reference genbank"""

    #     test_dir = os.path.join(GF_TEST_DIR, 'tenaillon',
    #             'Line20')
    #     variant_set = self._run_contig_walk_test(test_dir)

    #     me_variants = []
    #     for variant in variant_set.variants.all():
    #         for vccd in variant.variantcallercommondata_set.all():
    #             if vccd.data.get('INFO_METHOD', None) == 'ME_GRAPH_WALK':
    #                 me_variants.append(variant)

    #     self.assertTrue(me_variants)
    #     for variant in me_variants:
    #         self.assertTrue(1305000 < variant.position < 1306000)
    #         alts = variant.get_alternates()
    #         self.assertTrue(len(alts) == 1)
    #         alt = alts[0]
    #         self.assertTrue(1000 < len(alt) < 2000)
