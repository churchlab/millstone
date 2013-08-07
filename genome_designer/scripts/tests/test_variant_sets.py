"""
Tests for adding and removing variants from variant_sets.
"""

from django.core.files.uploadedfile import UploadedFile
from django.test import TestCase

from main.models import Project
from main.models import ReferenceGenome
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Variant
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.variant_sets import add_or_remove_variants_from_set
from scripts.bootstrap_data import bootstrap_fake_data
from settings import PWD as GD_ROOT_PATH

from scripts.bootstrap_data import TEST_PROJECT_NAME
from scripts.bootstrap_data import REF_GENOME_1_LABEL
from scripts.bootstrap_data import SAMPLE_1_LABEL
from scripts.bootstrap_data import VARIANTSET_1_LABEL, VARIANTSET_2_LABEL

class TestAddAndRemoveVariantsFromSet(TestCase):
    """Tests for scripts.import_util.import_samples_from_targets_file().
    """

    def setUp(self):
        bootstrap_fake_data()

        self.ref_genome_1 = ReferenceGenome.objects.get(
            label=REF_GENOME_1_LABEL)

        (self.sample_1, created) = ExperimentSample.objects.get_or_create(
                project__title=TEST_PROJECT_NAME,
                label=SAMPLE_1_LABEL)

        self.var_set1_uid = VariantSet.objects.get_or_create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_1_LABEL)[0].uid

        self.var_set2_uid = VariantSet.objects.get_or_create(
                reference_genome=self.ref_genome_1,
                label=VARIANTSET_2_LABEL)[0].uid

    def test_add_variants_to_set(self):

        variant_uids = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                chromosome='chrom',
                position__gt=70).values_list('uid')

        response = add_or_remove_variants_from_set(
                variant_uids,
                'add',
                self.var_set1_uid)

        assert response['alert_type'] == 'info', str(response)

    def test_remove_variant_from_set(self):

        variant_uids = Variant.objects.filter(
                reference_genome=self.ref_genome_1,
                chromosome='chrom',
                variantset__uid=self.var_set2_uid).values_list('uid')

        print self.var_set1_uid
        print variant_uids

        response = add_or_remove_variants_from_set(
                variant_uids,
                'remove',
                self.var_set2_uid)

        assert response['alert_type'] == 'info', str(response)

