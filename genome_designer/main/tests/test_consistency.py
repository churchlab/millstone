
from django.test import TestCase


from main.consistency import ensure_variant_set_consistency
from main.models import AlignmentGroup
from main.models import ExperimentSample
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities


SAMPLE_1_LABEL = 'sample1'
SAMPLE_2_LABEL = 'sample2'
VARIANTSET_1_LABEL = 'New Set A'


class TestEnsureVariantSetConsistency(TestCase):

    def test_simple(self):

        common_entities = create_common_entities()
        project = common_entities['project']
        self.ref_genome_1 = common_entities['reference_genome']
        self.chromosome = common_entities['chromosome']

        self.sample_1 = ExperimentSample.objects.create(
                project=project,
                label=SAMPLE_1_LABEL)

        self.sample_2 = ExperimentSample.objects.create(
                project=project,
                label=SAMPLE_2_LABEL)

        alignment_group = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=self.ref_genome_1,
            aligner=AlignmentGroup.ALIGNER.BWA)

        var_set_1 = VariantSet.objects.create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_1_LABEL)

        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome_1,
                chromosome=self.chromosome,
                position=100,
                ref_value='A')
        variant.variantalternate_set.add(
                VariantAlternate.objects.create(
                        variant=variant,
                        alt_value='G'))

        vtvs = VariantToVariantSet.objects.create(
                    variant=variant,
                    variant_set=var_set_1)

        vccd = VariantCallerCommonData.objects.create(
                variant=variant,
                source_dataset_id=1,
                alignment_group=alignment_group,
                data={})

        # Add a VariantEvidence with no GT_TYPE.
        VariantEvidence.objects.create(
                experiment_sample=self.sample_1,
                variant_caller_common_data=vccd,
                data={})
        self.assertEqual(0, vtvs.sample_variant_set_association.count())

        # Now add a VariantEvidence that has GT_TYPE=2 and run
        # ensure_variant_set_consistency().
        raw_sample_data_dict = {
                'CALLED': True,
                'GT_TYPE': 2,
                'GT_BASES': 'G/G'
        }
        VariantEvidence.objects.create(
                experiment_sample=self.sample_2,
                variant_caller_common_data=vccd,
                data=raw_sample_data_dict)
        ensure_variant_set_consistency(var_set_1)
        self.assertEqual(1, vtvs.sample_variant_set_association.count())
