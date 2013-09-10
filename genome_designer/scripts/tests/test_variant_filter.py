"""
Tests for variant_filter.py.
"""

from django.contrib.auth.models import User
from django.test import TestCase

from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant

from scripts.variant_filter import get_variants_that_pass_filter

class TestVariantFilter(TestCase):

    def setUp(self):
        user = User.objects.create_user('testuser', password='password',
                email='test@test.com')
        project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=project,
                label='refgenome', num_chromosomes=1, num_bases=1000)


    def test_filter__by_position(self):
        """Test filtering by position.
        """
        # Create several Variants with positions:
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        for pos in range(10):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome='chrom',
                position=pos,
                ref_value='A',
                alt_value='G')

        # Test querying Variants with position > 5.
        variants_above_5 = get_variants_that_pass_filter('position > 5',
                self.ref_genome)
        self.assertEqual(4, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var.position > 5)

        # Test querying Variants with position >= 5.
        variants_above_5 = get_variants_that_pass_filter('position >= 5',
                self.ref_genome)
        self.assertEqual(5, len(variants_above_5))
        for var in variants_above_5:
            self.assertTrue(var.position >= 5)


    def test_filter__by_chromosome(self):
        """Test filtering by position.
        """
        CHROM_1 = 'chrom'
        for pos in range(6):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_1,
                position=pos,
                ref_value='A',
                alt_value='G')

        CHROM_2 = 'chrom2'
        for pos in range(9):
            Variant.objects.create(
                type=Variant.TYPE.TRANSITION,
                reference_genome=self.ref_genome,
                chromosome=CHROM_2,
                position=pos,
                ref_value='A',
                alt_value='G')

        self.assertEqual(6, len(get_variants_that_pass_filter(
                'chromosome = chrom', self.ref_genome)))

        self.assertEqual(9, len(get_variants_that_pass_filter(
                'chromosome = chrom2', self.ref_genome)))
