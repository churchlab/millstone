"""
Tests for Django views.
"""

from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.test import Client
from django.test import TestCase

from main.models import Project
from scripts.bootstrap_data import bootstrap_fake_data
from scripts.bootstrap_data import TEST_USERNAME
from scripts.bootstrap_data import TEST_PASSWORD


STATUS_CODE__SUCCESS = 200
STATUS_CODE__NOT_LOGGED_IN = 302
STATUS_CODE__NOT_FOUND = 404


class TestViews(TestCase):

    def setUp(self):
        bootstrap_fake_data()

        test_project = Project.objects.all()[0]
        ref_genome = test_project.referencegenome_set.get(label='test_genome')
        alignment_group = ref_genome.alignmentgroup_set.all()[0]

        # Urls that do not require the user to be logged in.
        self.no_login_required_urls = [
                reverse('genome_designer.main.views.home_view'),
        ]

        # Urls that require the user to be logged in, but do not try any
        # particular entity.
        self.non_specific_login_required_urls = [
                reverse('genome_designer.main.views.project_list_view'),
                reverse('genome_designer.main.views.project_create_view'),
        ]

        # Urls for a specific entity.
        self.specific_entity_urls = [
                # Project-specific views
                reverse('genome_designer.main.views.project_view',
                        args=(test_project.uid,)),

                # Reference genomes
                reverse('genome_designer.main.views.reference_genome_list_view',
                        args=(test_project.uid,)),
                reverse('genome_designer.main.views.reference_genome_view',
                        args=(test_project.uid, ref_genome.uid)),

                # Alignments
                reverse('genome_designer.main.views.alignment_list_view',
                        args=(test_project.uid,)),
                reverse('genome_designer.main.views.alignment_create_view',
                        args=(test_project.uid,)),
                reverse('genome_designer.main.views.alignment_view',
                        args=(test_project.uid, alignment_group.uid)),

                # Variant sets
                reverse('genome_designer.main.views.variant_set_list_view',
                        args=(test_project.uid,)),

                # Samples
                reverse('genome_designer.main.views.sample_list_view',
                        args=(test_project.uid,)),

                # Variants
                reverse('genome_designer.main.views.variant_list_view',
                        args=(test_project.uid,)),

                # Genes
                reverse('genome_designer.main.views.gene_list_view',
                        args=(test_project.uid,)),

                # GO terms
                reverse('genome_designer.main.views.goterm_list_view',
                        args=(test_project.uid,)),
        ]

        # The fake web browser client used to make requests.
        self.client = Client()

    def assert_url_response(self, url, expected_status_code):
        """Helper method that calls a URL and compares the response status
        code to expected_status_code.
        """
        response = self.client.get(url)
        self.assertEqual(expected_status_code, response.status_code,
                ("Simple url test failed for %s with status code %d. " +
                        "Expected status code %d.") % (
                                url, response.status_code, expected_status_code))


    def test_views__logged_out(self):
        """Tests calling the views without a logged in user.
        """
        login_error_urls = (self.non_specific_login_required_urls +
                self.specific_entity_urls)
        for url in login_error_urls:
            self.assert_url_response(url, STATUS_CODE__NOT_LOGGED_IN)

        success_urls = self.no_login_required_urls
        for url in success_urls:
            self.assert_url_response(url, STATUS_CODE__SUCCESS)


    def test_views__logged_in_owner(self):
        """Tests calling views with the owner logged in.
        """
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)
        all_urls = (self.no_login_required_urls +
                self.non_specific_login_required_urls +
                self.specific_entity_urls)
        for url in all_urls:
            self.assert_url_response(url, STATUS_CODE__SUCCESS)


    def test_views__logged_in_non_owner(self):
         """Tests calling views with the non-owner logged in.
         """
         OTHER_USERNAME = 'justtest'
         OTHER_PASSWORD = 'other_password'
         OTHER_EMAIL = 'justtest@me.com'
         user = User.objects.create_user(
                 OTHER_USERNAME, password=OTHER_PASSWORD, email=OTHER_EMAIL)
         self.client.login(username=OTHER_USERNAME, password=OTHER_PASSWORD)

         error_urls = self.specific_entity_urls
         for url in error_urls:
             self.assert_url_response(url, STATUS_CODE__NOT_FOUND)

         success_urls = (self.non_specific_login_required_urls +
                 self.no_login_required_urls)
         for url in success_urls:
             self.assert_url_response(url, STATUS_CODE__SUCCESS)