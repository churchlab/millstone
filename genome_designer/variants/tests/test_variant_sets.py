"""
Tests for adding and removing variants from variant_sets.
"""

import os
import random

from django.conf import settings
from django.test import TestCase
from interval import interval

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantSet
from main.testing_util import create_common_entities
from scripts.bootstrap_data import create_fake_variants_and_variant_sets
from scripts.import_util import add_dataset_to_entity
from scripts.import_util import copy_dataset_to_entity_data_dir
from variants.variant_sets import add_variants_to_set_from_bed
from variants.variant_sets import MODIFY_VARIANT_SET_MEMBERSHIP__ADD
from variants.variant_sets import MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE
from variants.variant_sets import update_variant_in_set_memberships


TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')
TEST_BED = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'bed_test.bed')
SAMPLE_1_LABEL = 'sample1'
VARIANTSET_1_LABEL = 'New Set A'
VARIANTSET_2_LABEL = 'New Set B'


class TestAddVariantsToSetFromBed(TestCase):

    def test_add_variants_to_set_from_bed(self):

        common_entities = create_common_entities()
        project = common_entities['project']
        self.ref_genome_1 = common_entities['reference_genome']

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.ref_genome_1,
            aligner=AlignmentGroup.ALIGNER.BWA)

        (self.sample_1, created) = ExperimentSample.objects.get_or_create(
                project=project,
                label=SAMPLE_1_LABEL)

        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=self.sample_1)

        # Create variants in the bed regions from best_test.bed
        for var_poor_map in range(20):
            variant = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome_1,
                    chromosome='Chromosome',
                    position=random.randint(101,200),
                    ref_value='A')

            vccd = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset_id=1,
                alignment_group=alignment_group,
                data={}
            )

        for var_no_cov in range(20):
            variant = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome_1,
                    chromosome='Chromosome',
                    position=random.randint(301,400),
                    ref_value='A')

            vccd = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset_id=1,
                alignment_group=alignment_group,
                data={}
            )

            variant = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome_1,
                    chromosome='Chromosome',
                    position=random.randint(501,600),
                    ref_value='A')

            vccd = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset_id=1,
                alignment_group=alignment_group,
                data={}
            )

        new_bed_path = copy_dataset_to_entity_data_dir(
                entity= sample_alignment, 
                original_source_location= TEST_BED)

        bed_dataset = add_dataset_to_entity(sample_alignment, 
                dataset_label= Dataset.TYPE.BED_CALLABLE_LOCI,
                dataset_type= Dataset.TYPE.BED_CALLABLE_LOCI,
                filesystem_location= new_bed_path)

        vs_to_v_map = add_variants_to_set_from_bed(
                sample_alignment, bed_dataset)

        for variant_set, variants in vs_to_v_map.items():
            for v in variants:
                # POOR MAPPING QUAL should be from 101 to 200
                if variant_set.label == 'POOR_MAPPING_QUALITY':
                    self.assertTrue(v.position in interval(
                            [101, 200]))
                # NO COVERAGE should be from 301 to 400, 501 to 600
                elif variant_set.label == 'NO_COVERAGE':
                    self.assertTrue(v.position in interval(
                            [301,400],[501,600]))
                else:
                    raise AssertionError(
                            'bad variant set %s made.' % variant_set.label)

class TestAddAndRemoveVariantsFromSet(TestCase):

    def setUp(self):
        common_entities = create_common_entities()
        project = common_entities['project']
        self.ref_genome_1 = common_entities['reference_genome']

        create_fake_variants_and_variant_sets(self.ref_genome_1)

        (self.sample_1, created) = ExperimentSample.objects.get_or_create(
                project=project,
                label=SAMPLE_1_LABEL)

        self.var_set1 = VariantSet.objects.create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_1_LABEL)

        self.var_set2 = VariantSet.objects.create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_2_LABEL)

    def test_add_and_remove(self):
        """Test add and remove without samples involved.

        TODO: Make generating Variants deterministic. Right now we use
        random positions.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        variant_obj_list = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__gt=25,
                chromosome='chrom')
        variant_uid_sample_uid_pair_list = [obj.uid
                for obj in variant_obj_list]

        ### Test add.
        response = update_variant_in_set_memberships(
                self.ref_genome_1,
                variant_uid_sample_uid_pair_list,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid)
        self.assertEqual(response['alert_type'], 'info', str(response))
        self.assertEqual(len(variant_obj_list),
                self.var_set1.variants.all().count())

        ### Test remove.
        num_variants_in_set1_before_remove = (
                self.var_set1.variants.all().count())
        self.assertTrue(self.var_set1.variants.all().count() > 0)
        variant_obj_list = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__gt=75,
                chromosome='chrom',
                variantset__uid=self.var_set1.uid)
        variant_uids_to_remove = [obj.uid for obj in variant_obj_list]
        self.assertTrue(len(variant_uids_to_remove) > 0)

        response = update_variant_in_set_memberships(
                self.ref_genome_1,
                variant_uids_to_remove,
                MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE,
                self.var_set1.uid)
        self.assertEqual(response['alert_type'], 'info', str(response))

        # Check that we've reduced the number of variants.
        self.assertEqual(
                num_variants_in_set1_before_remove -
                        len(variant_uids_to_remove),
                self.var_set1.variants.all().count())
