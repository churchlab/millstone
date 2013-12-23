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
from main.models import VariantAlternate
from main.models import VariantCallerCommonData
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from scripts.vcf_parser import parse_alignment_group_vcf
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
        num_experiment_samples = len(experiment_sample_uids)
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


        variant_list = Variant.objects.filter(
                reference_genome=self.reference_genome)

        # There should be one Variant object for each record.
        self.assertEqual(record_count, len(variant_list))

        # Spot-check a few variants.
        self.assertEqual(1, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=376)))

        v_453 = Variant.objects.get(reference_genome=self.reference_genome,
                position=453)
        self.assertEqual(['G'], v_453.get_alternates())

        # Check false negatives.
        self.assertEqual(0, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=454)))

        # There should be one VariantCallerCommonData object for each record.
        self.assertEqual(record_count,
                len(VariantCallerCommonData.objects.filter(
                        variant__reference_genome=self.reference_genome)))

        # There should also be one VariantEvidence object per Variant x Sample.
        for variant in variant_list:
            vccd = variant.variantcallercommondata_set.all()[0]
            self.assertEqual(num_experiment_samples,
                    len(vccd.variantevidence_set.all()))

        # Check that alternate data is populated.
        #Chromosome  1330    .   CG  C,GC,AG 126.036 .   AB=0.5,0.5,1;ABP=3.0103,3.0103,7.35324;AC=1,1,1;AF=0.0833333,0.0833333,0.0833333;AN=12;AO=1,1,2;CIGAR=1M1D,2X,1X1M;DP=10;DPRA=1.33333,1.33333,1.33333;EPP=5.18177,5.18177,3.0103;EPPR=4.45795;HWE=-16.5861;LEN=1,2,1;MEANALT=2,2,1;MQM=60,37,48.5;MQMR=40.8333;NS=6;NUMALT=3;ODDS=1.50408;PAIRED=1,0,0.5;PAIREDR=0.166667;RO=6;RPP=5.18177,5.18177,7.35324;RPPR=16.0391;RUN=1,1,1;SAP=5.18177,5.18177,3.0103;SRP=4.45795;TYPE=del,mnp,snp;XAI=0,0.0102041,0.00515464;XAM=0,0.0102041,0.0253649;XAS=0,0,0.0202103;XRI=0.0016835;XRM=0.00835084;XRS=0.00666733;technology.illumina=1,1,1;BVAR GT:DP:RO:QR:AO:QA:GL    .   0/0:1:1:36:0,0,0:0,0,0:0,-0.30103,-3.6,-0.30103,-3.6,-3.6,-0.30103,-3.6,-3.6,-3.6   0/0:2:2:76:0,0,0:0,0,0:0,-0.60206,-7.03,-0.60206,-7.03,-7.03,-0.60206,-7.03,-7.03,-7.03 1/2:2:0:0:1,1,0:108,31,0:-8.645,-3.40103,-3.1,-6.30103,-0.30103,-6,-8.645,-3.40103,-6.30103,-8.645  .   0/3:2:0:0:0,0,2:0,0,73:-6.935,-6.935,-6.935,-6.935,-6.935,-6.935,-0.60206,-0.60206,-0.60206,0   0/0:2:2:72:0,0,0:0,0,0:0,-0.60206,-6.84,-0.60206,-6.84,-6.84,-0.60206,-6.84,-6.84,-6.84 .   0/0:1:1:34:0,0,0:0,0,0:0,-0.30103,-3.4,-0.30103,-3.4,-3.4,-0.30103,-3.4,-3.4,-3.4   .
        v_1330 = Variant.objects.get(reference_genome=self.reference_genome,
                position=1330)
        self.assertEqual(set(v_1330.get_alternates()),set(['C','GC','AG']))
        v_1330_c = VariantAlternate.objects.get(variant=v_1330, alt_value='C')
        self.assertTrue(len(v_1330_c.variantevidence_set.all()))
        v_1330_gc = VariantAlternate.objects.get(variant=v_1330, alt_value='GC')
        self.assertTrue(len(v_1330_gc.variantevidence_set.all()))
        self.assertEqual(v_1330_c.as_dict()['INFO_ABP'],
                v_1330_gc.as_dict()['INFO_ABP'])
