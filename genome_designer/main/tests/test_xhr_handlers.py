"""
Tests for xhr_handlers.py.
"""

import json
import os
import pickle

from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.http.request import HttpRequest
from django.test import Client
from django.test import TestCase

from main.adapters import OBJ_LIST
from main.models import Dataset
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantCallerCommonData
from main.xhr_handlers import _determine_visible_field_names
from main.xhr_handlers import VARIANT_LIST_REQUEST_KEY__PROJECT_UID
from main.xhr_handlers import VARIANT_LIST_REQUEST_KEY__REF_GENOME_UID
from main.xhr_handlers import VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS
from main.xhr_handlers import VARIANT_LIST_RESPONSE_KEY__LIST
from main.xhr_handlers import VARIANT_LIST_RESPONSE_KEY__TOTAL
from main.xhr_handlers import VARIANT_LIST_RESPONSE_KEY__SET_LIST
from main.xhr_handlers import VARIANT_LIST_RESPONSE_KEY__KEY_MAP
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.dynamic_snp_filter_key_map import update_filter_key_map
from settings import PWD as GD_ROOT


TEST_USERNAME = 'testusername'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@test.com'

TEST_DIR = os.path.join(GD_ROOT, 'test_data', 'genbank_aligned')
TEST_ANNOTATED_VCF = os.path.join(TEST_DIR, 'bwa_align_annotated.vcf')

STATUS_CODE__NOT_FOUND = 404
STATUS_CODE__NOT_LOGGED_IN = 302
STATUS_CODE__SUCCESS = 200

class TestGetVariantList(TestCase):

    url = reverse('genome_designer.main.xhr_handlers.get_variant_list')

    def setUp(self):
        # Useful models.
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=self.project,
                label='refgenome', num_chromosomes=1, num_bases=1000)

        # Make sure the reference genome has the required vcf keys.
        initialize_filter_key_map(self.ref_genome)
        update_filter_key_map(self.ref_genome, TEST_ANNOTATED_VCF)
        self.vcf_dataset = Dataset.objects.create(
                label='test_data_set',
                type=Dataset.TYPE.VCF_FREEBAYES,
                filesystem_location=TEST_ANNOTATED_VCF)


        # Fake web browser client used to make requests.
        self.client = Client()
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)


    def test__logged_out(self):
        """Test that logged out fails.
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertEqual(STATUS_CODE__NOT_LOGGED_IN, response.status_code)


    def test__missing_params(self):
        response = self.client.get(self.url)
        self.assertEqual(STATUS_CODE__NOT_FOUND, response.status_code)


    def test__basic_function(self):
        """Basic test.
        """
        TOTAL_NUM_VARIANTS = 10
        for pos in range(TOTAL_NUM_VARIANTS):
            Variant.objects.create(
                    type=Variant.TYPE.TRANSITION,
                    reference_genome=self.ref_genome,
                    chromosome='chrom',
                    position=pos,
                    ref_value='A')

        request_data = {
            'refGenomeUid': self.ref_genome.uid,
            'projectUid': self.project.uid
        }
        response = self.client.get(self.url, request_data)

        self.assertEqual(STATUS_CODE__SUCCESS, response.status_code)

        EXPECTED_RESPONSE_KEYS = [
            VARIANT_LIST_RESPONSE_KEY__LIST,
            VARIANT_LIST_RESPONSE_KEY__TOTAL,
            VARIANT_LIST_RESPONSE_KEY__SET_LIST,
            VARIANT_LIST_RESPONSE_KEY__KEY_MAP,
        ]

        response_data = json.loads(response.content)

        # Make sure expected keys in response.
        if not all(key in response_data for key in EXPECTED_RESPONSE_KEYS):
            self.fail('Missing keys in response')

        # TODO(gleb): Uncomment when we fix total counts.
        # self.assertEqual(TOTAL_NUM_VARIANTS,
        #         response_data[VARIANT_LIST_RESPONSE_KEY__TOTAL])


    def test_determine_visible_field_names(self):
        EXPECTED_VISIBLE_KEYS = ['INFO_EFF_EFFECT']
        request = HttpRequest()
        request.GET[VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS] = json.dumps(
                EXPECTED_VISIBLE_KEYS)
        self.assertEqual(EXPECTED_VISIBLE_KEYS,
                _determine_visible_field_names(request, '', self.ref_genome))


    # def test__INFO_EFF(self):
    #     """Basic test.
    #     """
    #     variant = Variant.objects.create(
    #             type=Variant.TYPE.TRANSITION,
    #             reference_genome=self.ref_genome,
    #             chromosome='chrom',
    #             position=100,
    #             ref_value='A')

    #     raw_data_dict = {
    #             'INFO_EFF_EFFECT': pickle.dumps(['NON_SYNONYMOUS'])
    #     }
    #     VariantCallerCommonData.objects.create(
    #             variant=variant,
    #             source_dataset=self.vcf_dataset,
    #             data=raw_data_dict
    #     )

    #     VISIBLE_KEYS = ['INFO_EFF_EFFECT']

    #     request_data = {
    #         VARIANT_LIST_REQUEST_KEY__REF_GENOME_UID: self.ref_genome.uid,
    #         VARIANT_LIST_REQUEST_KEY__PROJECT_UID: self.project.uid,
    #         VARIANT_LIST_REQUEST_KEY__VISIBLE_KEYS: json.dumps(VISIBLE_KEYS)
    #     }

    #     response = self.client.get(self.url, request_data)

    #     self.assertEqual(STATUS_CODE__SUCCESS, response.status_code)

    #     response_data = json.loads(response.content)

    #     # TODO(gleb): Uncomment when we fix total counts.
    #     # self.assertEqual(1,
    #     #         response_data[VARIANT_LIST_RESPONSE_KEY__TOTAL])

    #     # NOTE: This double layer of json is not good.
    #     fe_variant = json.loads(
    #             response_data[VARIANT_LIST_RESPONSE_KEY__LIST])[OBJ_LIST][0]

    #     # NOTE: This is the current output. If you break this test and are
    #     # doing something better, please change the test to account for the
    #     # new behavior.
    #     self.assertEqual("['NON_SYNONYMOUS']", fe_variant['INFO_EFF_EFFECT'])
