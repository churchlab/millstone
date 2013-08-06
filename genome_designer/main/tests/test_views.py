"""
Tests for Django views.
"""

from django.core.urlresolvers import reverse
from django.test import Client
from django.test import TestCase

from main.models import Project
from scripts.bootstrap_data import bootstrap_fake_data
from scripts.bootstrap_data import TEST_USERNAME
from scripts.bootstrap_data import TEST_PASSWORD

STATIC_CODE__SUCCESS = 200

class TestViews(TestCase):

    def setUp(self):
        bootstrap_fake_data()

        # The fake web browser client used to make requests.
        self.client = Client()

        # Login for all tests for now.
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)


    def test_views_simple(self):
        """Call all of the views and make sure they don't error.
        """
        test_project = Project.objects.all()[0]
        ref_genome = test_project.referencegenome_set.get(label='test_genome')
        alignment_group = ref_genome.alignmentgroup_set.all()[0]

        urls = [
                reverse('genome_designer.main.views.home_view'),

                # Project-specific views
                reverse('genome_designer.main.views.project_list_view'),
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
        for url in urls:
            response = self.client.get(url)
            self.assertEqual(STATIC_CODE__SUCCESS,
                    response.status_code,
                    "Simple url test failed for %s with status code %d" % (
                        url, response.status_code))
