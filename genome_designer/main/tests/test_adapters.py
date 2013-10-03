"""
Tests for adapters.py
"""

import json

from django.contrib.auth.models import User
from django.test import TestCase

from main.adapters import adapt_model_to_frontend
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant

class TestAdapters(TestCase):

    def setUp(self):
        """Override.
        """
        TEST_USERNAME = 'testuser'
        TEST_PASSWORD = 'password'
        TEST_EMAIL = 'test@example.com'
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)

        TEST_PROJECT_NAME = 'recoli'
        test_project = Project.objects.create(
            title=TEST_PROJECT_NAME,
            owner=user.get_profile())

        REF_GENOME_1_LABEL = 'mg1655'
        self.ref_genome_1 = ReferenceGenome.objects.create(
            label=REF_GENOME_1_LABEL, project=test_project, num_chromosomes=1,
            num_bases=100)


    def test_adapters__one_level(self):
        """Test adapting to a single level.
        """
        fe_ref_genomes = json.loads(adapt_model_to_frontend(
                ReferenceGenome, {'id': self.ref_genome_1.id}))

        self.assertTrue('field_config' in fe_ref_genomes)
        self.assertTrue('obj_list' in fe_ref_genomes)
        self.assertEqual(1, len(fe_ref_genomes['obj_list']))
        ref_genome_1_fe = fe_ref_genomes['obj_list'][0]
        for field in ReferenceGenome.get_field_order():
            self.assertTrue(field['field'] in ref_genome_1_fe)
        self.assertTrue('href' in ref_genome_1_fe)


    def test_adapters__two_levels(self):
        """Test adapting model with nested models.
        """
        variant = Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome_1,
                chromosome='chrom',
                position=100,
                ref_value='A',
                alt_value='G')

        fe_variants = json.loads(adapt_model_to_frontend(Variant))

        self.assertEqual(1, len(fe_variants['obj_list']))
        fe_variant = fe_variants['obj_list'][0]

        for field in Variant.get_field_order():
            self.assertTrue(field['field'] in fe_variant)

        # self.assertTrue('reference_genome' in fe_variant)
        
        # for ref_genome_field in ReferenceGenome.get_field_order():
        #     self.assertTrue(ref_genome_field['field'] in 
        #         fe_variant['reference_genome'])
