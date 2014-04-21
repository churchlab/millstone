import csv
import StringIO

from django.core.urlresolvers import reverse
from django.test import Client
from django.test import TestCase

from main.upload_template_views import SAMPLE_LIST_BROWSER_UPLOAD_TEMPLATE
from utils.import_util import SAMPLE_BROWSER_UPLOAD_KEY__SAMPLE_NAME
from utils.import_util import SAMPLE_BROWSER_UPLOAD_KEY__READ_1
from utils.import_util import SAMPLE_BROWSER_UPLOAD_KEY__READ_2


class TestSampleListBrowserUploadTemplate(TestCase):

    def setUp(self):
        self.client = Client()

    def test_template(self):
        """Makes sure the template is as expected.

        The main issue we're testing here is that we don't accidentally
        break our template.
        """
        # Make the request.
        template_url  = reverse(
                'main.upload_template_views.sample_list_browser_upload_template')
        response = self.client.get(template_url)

        # Create a reader on the response content.
        reader = csv.DictReader(StringIO.StringIO(response.content),
                delimiter='\t')

        # Check that the expected headrs are there.
        targets_file_header = reader.fieldnames
        EXPECTED_TEMPLATE_HEADER = [
                SAMPLE_BROWSER_UPLOAD_KEY__SAMPLE_NAME,
                SAMPLE_BROWSER_UPLOAD_KEY__READ_1,
                SAMPLE_BROWSER_UPLOAD_KEY__READ_2
        ]

        missing_header_cols = (set(EXPECTED_TEMPLATE_HEADER) -
            set(targets_file_header))
        self.assertEqual(0, len(missing_header_cols),
            "Missing cols: %s. "
            "Make sure the template is properly tab-delimited." %
            ' '.join(missing_header_cols))
