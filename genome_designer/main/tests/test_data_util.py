"""
Tests for data_util.py
"""

from django.contrib.auth.models import User
from django.test import TestCase


from main.data_util import RequestScopedVariantDataCache
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant


class TestRequestScopedVariantDataCache(TestCase):

    def setUp(self):
        user = User.objects.create_user('testuser', password='password',
                email='test@test.com')
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=self.project,
                label='refgenome', num_chromosomes=1, num_bases=1000)


    def test_populate(self):
        """Tests populate.

        Also serves as means of profiling requests.
        """
        from django.conf import settings

        variant_list = []
        for pos in range(10):
            variant_list.append(Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A'))

        # Just testing constructor and populate for now.
        # TODO: Add more thorough tests if we stick with this design.
        data_cache = RequestScopedVariantDataCache()
        data_cache.populate(variant_list, self.project)
