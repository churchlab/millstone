"""Tests for reference_genome_maker_util.py.
"""

import os

from Bio import SeqIO
from django.conf import settings
from django.test import TestCase

from main.model_utils import get_dataset_with_type
from main.models import Chromosome
from main.models import Dataset
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities
from utils.import_util import import_reference_genome_from_local_file
from utils.reference_genome_maker_util import generate_new_reference_genome


TEST_DATA_DIR = os.path.join(settings.PWD, 'test_data')
TEST_GENBANK = os.path.join(TEST_DATA_DIR, 'full_vcf_test_set',
        'mg1655_tolC_through_zupT.gb')


class TestReferenceGenomeMakerUtil(TestCase):

    def setUp(self):
        """Override.
        """
        self.common_entities = create_common_entities()
        self.project = self.common_entities['project']

    def test_basic(self):
        """Basic test.
        """
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_GENBANK, 'genbank')
        variant_set = VariantSet.objects.create(
                reference_genome=self.reference_genome,
                label='vs1')

        ref_genome_filepath = get_dataset_with_type(self.reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

        with open(ref_genome_filepath) as fh:
            ref_genome_seq_record = SeqIO.read(fh, 'genbank')

        for position in range(10, 111, 10):
            ref_value = ref_genome_seq_record[position - 1]
            var = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.reference_genome,
                    chromosome=Chromosome.objects.get(reference_genome=self.reference_genome),
                    position=position,
                    ref_value=ref_value)

            VariantAlternate.objects.create(
                    variant=var, alt_value='G')

            VariantToVariantSet.objects.create(
                    variant=var, variant_set=variant_set)

        new_ref_genome_params = {
            'label': 'new'
        }

        new_ref_genome = generate_new_reference_genome(
                variant_set, new_ref_genome_params)

        new_ref_genome_filepath = get_dataset_with_type(
                        new_ref_genome,
                        Dataset.TYPE.REFERENCE_GENOME_GENBANK)\
                .get_absolute_location()
        with open(new_ref_genome_filepath) as fh:
            new_ref_genome_seq_record = SeqIO.read(fh, 'genbank')

        # Assert size unchangd.
        self.assertEqual(len(new_ref_genome_seq_record),
                len(ref_genome_seq_record))

        # Assert mutations are there.
        for position in range(10, 111, 10):
            self.assertEqual('G', str(new_ref_genome_seq_record[position - 1]))

        # Assert new genome is annotated.
        new_ref_genome.is_annotated()