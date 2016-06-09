import os

from django.conf import settings
from django.test import TestCase

from main.model_utils import get_dataset_with_type
from main.models import Chromosome
from main.models import Dataset
from main.testing_util import create_common_entities
from utils.import_util import import_reference_genome_from_local_file


TEST_DATA_DIR = os.path.join(settings.PWD, 'test_data')
TEST_GENBANK = os.path.join(TEST_DATA_DIR, 'mg1655.genbank')


class TestGenbankUtil(TestCase):

    def setUp(self):
        """Override.
        """
        self.common_entities = create_common_entities()
        self.project = self.common_entities['project']

    def test_generate_genbank_mobile_element_multifasta(self):
        """Test generation of the mobile element fasta.
        """
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_GENBANK, 'genbank')
        self.reference_genome.ensure_mobile_element_multifasta()

        me_fa_dataset = get_dataset_with_type(
                self.reference_genome,
                Dataset.TYPE.MOBILE_ELEMENT_FASTA)

        assert os.path.exists(
                me_fa_dataset.get_absolute_location())