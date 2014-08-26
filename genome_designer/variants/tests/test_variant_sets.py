"""
Tests for adding and removing variants from variant_sets.
"""

import os
import random

from django.conf import settings
from django.test import TestCase
import pyinter

from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities
from utils.import_util import add_dataset_to_entity
from utils.import_util import copy_dataset_to_entity_data_dir
from variants.variant_sets import add_variants_to_set_from_bed
from variants.variant_sets import MODIFY_VARIANT_SET_MEMBERSHIP__ADD
from variants.variant_sets import MODIFY_VARIANT_SET_MEMBERSHIP__REMOVE
from variants.variant_sets import update_variant_in_set_memberships
from variants.variant_sets import update_variant_in_set_memberships__all_matching_filter


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
                    chromosome=Chromosome.objects.get(reference_genome=self.ref_genome_1),
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
                    chromosome=Chromosome.objects.get(reference_genome=self.ref_genome_1),
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
                    chromosome=Chromosome.objects.get(reference_genome=self.ref_genome_1),
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

        variant_set_labels = set([vs.label for vs in vs_to_v_map.keys()])
        self.assertEqual(set(['POOR_MAPPING_QUALITY', 'NO_COVERAGE']),
                variant_set_labels)

        for variant_set, variants in vs_to_v_map.items():
            for v in variants:
                # POOR MAPPING QUAL should be from 101 to 200
                if variant_set.label == 'POOR_MAPPING_QUALITY':
                    self.assertTrue(v.position in pyinter.closedopen(
                            101, 200))
                # NO COVERAGE should be from 301 to 400, 501 to 600
                elif variant_set.label == 'NO_COVERAGE':
                    self.assertTrue(v.position in pyinter.IntervalSet([
                                    pyinter.closedopen(301,400),
                                    pyinter.closedopen(501,600)]))
                else:
                    raise AssertionError(
                            'bad variant set %s made.' % variant_set.label)


class TestAddAndRemoveVariantsFromSet(TestCase):

    def setUp(self):
        common_entities = create_common_entities()
        project = common_entities['project']
        self.ref_genome_1 = common_entities['reference_genome']
        self.chromosome = common_entities['chromosome']

        self.sample_1 = ExperimentSample.objects.create(
                project=project,
                label=SAMPLE_1_LABEL)

        # Create 100 variants with alts, positions 1 to 100.
        # We add them to a fake VariantSet so that they show up in the
        # materialized variant view table.
        self.first_variant_position = 1
        self.num_variants = 100
        fake_variant_set = VariantSet.objects.create(
                reference_genome=self.ref_genome_1,
                label='fake')

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.ref_genome_1,
            aligner=AlignmentGroup.ALIGNER.BWA)

        for position in xrange(self.first_variant_position,
                self.first_variant_position + self.num_variants):
            variant = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome_1,
                    chromosome=self.chromosome,
                    position=position,
                    ref_value='A')
            variant.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=variant,
                            alt_value='G'))
            vccd = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset_id=1,
                alignment_group=alignment_group,
                data={})

            VariantEvidence.objects.create(
                experiment_sample=self.sample_1,
                variant_caller_common_data=vccd,
                data={})

            VariantToVariantSet.objects.create(
                    variant=variant,
                    variant_set=fake_variant_set)

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
                chromosome=self.chromosome)
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
                chromosome=self.chromosome,
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

    def test_add__sample_association(self):
        """Tests adding a Variant to a VariantSet, with an association to a
        particular ExperimentSample.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        # Add these variants without sample associations.
        variants_to_add_with_no_sample_association = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__lt=25,
                chromosome=self.chromosome)
        variants_no_association_data_str_list = [obj.uid
                for obj in variants_to_add_with_no_sample_association]

        # Add these variants associated with sample 1.
        variants_to_add_with_sample_association = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                position__gte=25,
                chromosome=self.chromosome)
        variants_with_association_data_str_list = [
                obj.uid + ',' + self.sample_1.uid
                for obj in variants_to_add_with_sample_association]

        # Put these together and run it all the at the same time.
        all_data_str_list = (variants_no_association_data_str_list +
                variants_with_association_data_str_list)
        response = update_variant_in_set_memberships(
                self.ref_genome_1,
                all_data_str_list,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid)
        self.assertEqual(response['alert_type'], 'info', str(response))

        # Make sure all the variants are there.
        self.assertEqual(self.num_variants,
                self.var_set1.variants.all().count())

        all_vtvs = VariantToVariantSet.objects.filter(
                variant_set=self.var_set1)

        # Sanity check.
        self.assertEqual(self.num_variants, all_vtvs.count())

        # Check that the Variants we expected to have an association have it.
        for vtvs in all_vtvs:
            if vtvs.variant in variants_to_add_with_sample_association:
                self.assertEqual(1, vtvs.sample_variant_set_association.count())
                self.assertEqual(self.sample_1,
                        vtvs.sample_variant_set_association.all()[0])
            else:
                self.assertEqual(0, vtvs.sample_variant_set_association.count())

    def test_all_matching_filter__all__cast(self):
        """Test adding all matching '' filter, cast.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        # Add all variants.
        update_variant_in_set_memberships__all_matching_filter(
                self.ref_genome_1,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid,
                '',
                False)

        self.assertEqual(self.num_variants,
                self.var_set1.variants.all().count())

        # Make sure no ExperimentSample association.
        for vtvs in VariantToVariantSet.objects.filter(
                variant_set=self.var_set1):
            self.assertEqual(0, vtvs.sample_variant_set_association.count())

    def test_all_matching_filter__partial__cast(self):
        """Test adding all matching partial filter, cast.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        # Add all variants.
        update_variant_in_set_memberships__all_matching_filter(
                self.ref_genome_1,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid,
                'position <= 50',
                False)

        self.assertEqual(50, self.var_set1.variants.all().count())

        # Make sure no ExperimentSample association.
        for vtvs in VariantToVariantSet.objects.filter(
                variant_set=self.var_set1):
            self.assertEqual(0, vtvs.sample_variant_set_association.count())

    def test_all_matching_filter__all__melted(self):
        """Test adding all matching '' filter, melted.

        In the melted case, we expect the variants to be associated with all
        samples.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        # Add all variants.
        update_variant_in_set_memberships__all_matching_filter(
                self.ref_genome_1,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid,
                '',
                True)

        self.assertEqual(self.num_variants,
                self.var_set1.variants.all().count())

        # Make sure no ExperimentSample association.
        all_vtvs = VariantToVariantSet.objects.filter(
                variant_set=self.var_set1)

        # Sanity check.
        self.assertEqual(self.num_variants, all_vtvs.count())

        # Check that the Variants we expected to have an association have it.
        for vtvs in all_vtvs:
            self.assertEqual(1, vtvs.sample_variant_set_association.count())
            self.assertEqual(self.sample_1,
                    vtvs.sample_variant_set_association.all()[0])

    def test_all_matching_filter__partial__melted(self):
        """Test adding all matching partial filter, melted.

        In the melted case, we expect the variants to be associated with all
        samples.
        """
        # No variants before adding.
        self.assertEqual(0, self.var_set1.variants.all().count())

        EXPECTED_NUM_VARIANTS = 50

        # Add all variants.
        update_variant_in_set_memberships__all_matching_filter(
                self.ref_genome_1,
                MODIFY_VARIANT_SET_MEMBERSHIP__ADD,
                self.var_set1.uid,
                'position <= %d' % EXPECTED_NUM_VARIANTS,
                True)

        self.assertEqual(EXPECTED_NUM_VARIANTS,
                self.var_set1.variants.all().count())

        # Make sure no ExperimentSample association.
        all_vtvs = VariantToVariantSet.objects.filter(
                variant_set=self.var_set1)

        # Sanity check.
        self.assertEqual(EXPECTED_NUM_VARIANTS, all_vtvs.count())

        # Check that the Variants we expected to have an association have it.
        for vtvs in all_vtvs:
            self.assertEqual(1, vtvs.sample_variant_set_association.count())
            self.assertEqual(self.sample_1,
                    vtvs.sample_variant_set_association.all()[0])
