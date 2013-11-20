import tempfile

from django.test import TestCase
from scripts.util import ensure_line_lengths
from scripts.util import SNPEFF_MIN_LINE_LEN


class TestUtil(TestCase):
    def test_ensure_line_lengths(self):
        test_tmp = tempfile.NamedTemporaryFile(delete=False)
        print >> test_tmp, "test1"
        print >> test_tmp, ""
        print >> test_tmp, "test2              "
        test_tmp.close()
        ensure_line_lengths(test_tmp.name)
        with open(test_tmp.name, 'r') as check_test_fh:
            for line in check_test_fh:
                assert len(line) == SNPEFF_MIN_LINE_LEN + 1
