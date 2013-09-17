"""
Tests for vcf_parser.py
"""

import os

from django.contrib.auth.models import User
from django.test import TestCase
import vcf

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import Variant
from main.models import VariantCallerCommonData
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from scripts.vcf_parser import parse_alignment_group_vcf
from scripts.vcf_parser import populate_common_data_eff
from settings import PWD as GD_ROOT

TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_GENOME_SNPS = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome_snps.vcf')


class TestVCFParser(TestCase):

    def setUp(self):
        # Test models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')


    def test_parser(self):
        """Basic tests for the parser.
        """
        VCF_DATATYPE = Dataset.TYPE.VCF_FREEBAYES
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)
        copy_and_add_dataset_source(alignment_group, VCF_DATATYPE,
                VCF_DATATYPE, TEST_GENOME_SNPS)

        # Create experiment sample objects having UIDs that correspond to those
        # in the vcf file. This is a bit "fake" in that the actual pipeline we
        # will be generating the vcf file from the samples (see add_groups()
        # stage of pipeline.
        with open(TEST_GENOME_SNPS) as fh:
            reader = vcf.Reader(fh)
            experiment_sample_uids = reader.samples
        for sample_uid in experiment_sample_uids:
            ExperimentSample.objects.create(
                uid=sample_uid,
                project=self.project,
                label='fakename:' + sample_uid,
                group='Plate 1',
                well='A01',
                num_reads=100,
            )

        # Count the number of records in the vcf file for testing.
        record_count = 0
        with open(TEST_GENOME_SNPS) as fh:
            for record in vcf.Reader(fh):
                record_count += 1

        # Parse the vcf
        parse_alignment_group_vcf(alignment_group, VCF_DATATYPE)

        # There should be one VariantCallerCommonData object for each record.
        self.assertEqual(record_count,
                len(VariantCallerCommonData.objects.filter(
                        variant__reference_genome=self.reference_genome)))

        # There should also be one Variant object for each record.
        self.assertEqual(record_count, len(Variant.objects.filter(
                reference_genome=self.reference_genome)))

        # Spot-check a few variants.
        self.assertEqual(1, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=376)))
        self.assertEqual(1, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=453)))

        # Check false negatives.
        self.assertEqual(0, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=454)))

    def test_populate_common_data_eff(self):
        """ Test the regex on a few snpeff field examples.
        """
        data_dict = {}

        # single eff
        populate_common_data_eff(''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC',
            '||CODING|b3038|1|1)')), data_dict)

        # eff with errors
        populate_common_data_eff(''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC',
            '||CODING|b3038|1|1|WARN_TEST|ERROR_TEST)')), data_dict)

        # multi-eff
        data_dict = populate_common_data_eff(''.join((
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC',
            '||CODING|b3038|1|1|WARN_TEST|ERROR_TEST),',
            'NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aGg|M239T|386|ygiC',
            '||CODING|b3038|1|1|ERROR_TEST|WARN_TEST)')), data_dict)

        self.assertEqual(data_dict['INFO_EFF_CONTEXT'],['aTg/aCg','aTg/aGg'])
