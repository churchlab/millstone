import StringIO
from tempfile import mkstemp

from django.conf import settings
from django.core.urlresolvers import reverse
from django.test import Client
from django.test import TestCase

from utils.import_util import determine_template_delimiter
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
        EXPECTED_TEMPLATE_HEADER = [
                SAMPLE_BROWSER_UPLOAD_KEY__SAMPLE_NAME,
                SAMPLE_BROWSER_UPLOAD_KEY__READ_1,
                SAMPLE_BROWSER_UPLOAD_KEY__READ_2
        ]

        # Make the request.
        template_url = reverse(
                'main.upload_template_views.sample_list_browser_upload_template')
        response = self.client.get(template_url)

        # Save the content to temp file to match interface of
        # determine_template_delimiter().
        _, temp_file_location = mkstemp(dir=settings.TEMP_FILE_ROOT)
        with open(temp_file_location, 'w') as temp_fh:
            temp_fh.write(StringIO.StringIO(response.content).read())

        delim = determine_template_delimiter(
                temp_file_location, EXPECTED_TEMPLATE_HEADER)

        self.assertEqual(',', delim)
