"""Tests for optmage_util.py.
"""

import csv
import os
import random
from StringIO import StringIO

from Bio import SeqIO
from django.conf import settings
from django.test import TestCase
from optmage.oligo_designer import OLIGO_TARGET_REQUIRED_PARAMS

from main.models import Chromosome
from main.models import Dataset
from main.models import Variant
from main.models import VariantAlternate
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.testing_util import create_common_entities
from utils.import_util import import_reference_genome_from_local_file
from utils.optmage_util import print_mage_oligos
from utils.optmage_util import ReplicationOriginParams


TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')


class TestOptmageUtil(TestCase):

    def setUp(self):
        self.common_entities = create_common_entities()
        self.ref_genome = import_reference_genome_from_local_file(
                self.common_entities['project'], 'ref_genome', TEST_FASTA,
                'fasta')
        ref_genome_source = self.ref_genome.dataset_set.get(
                type=Dataset.TYPE.REFERENCE_GENOME_FASTA)\
                        .get_absolute_location()
        with open(ref_genome_source) as fh:
            self.ref_genome_seq_record = SeqIO.read(fh, 'fasta')

    def test_print_mage_oligos(self):
        var_set_1 = VariantSet.objects.create(
            reference_genome=self.ref_genome,
            label='vs1')
        POSITION_RANGE = range(100, 1001, 100)
        for position in POSITION_RANGE:
            ref_value = self.ref_genome_seq_record.seq[position]
            alt_value = random.choice('ACGT')
            while alt_value == ref_value:
                alt_value = random.choice('ACGT')
            variant = Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome,
                    chromosome=Chromosome.objects.get(reference_genome=self.ref_genome),
                    position=position,
                    ref_value=ref_value)
            VariantAlternate.objects.create(
                    variant=variant,
                    alt_value=alt_value)
            VariantToVariantSet.objects.create(
                    variant=variant,
                    variant_set=var_set_1)
        output = StringIO()
        target_id_prefix = 'o_'
        print_mage_oligos(var_set_1, output, target_id_prefix,
                ReplicationOriginParams.from_defaults())

        reader = csv.DictReader(StringIO(output.getvalue()))
        all_rows = [row for row in reader]
        self.assertEqual(10, len(all_rows))

        # Check at least some reqiured keys are there.
        for row in all_rows:
            for key in OLIGO_TARGET_REQUIRED_PARAMS:
                self.assertTrue('target_id' in row)
                self.assertEqual(90, int(row['oligo_size']))

        start_position_set = set([int(row['start']) for row in all_rows])
        self.assertEqual(set(POSITION_RANGE), start_position_set)
