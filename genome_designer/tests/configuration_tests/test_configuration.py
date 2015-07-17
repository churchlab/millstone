"""
Tests configuration, e.g. proper paths, etc.
"""

import os

from django.conf import settings
from django.test import TestCase

from millstone_setup import download_tools
from conf import external_tool_registry


class TestConfiguration(TestCase):

    def test_tool_download(self):
        """Runs millstone_setup.py and makes sure all binaries are downloaded.

        NOTE: In development. Currently only tests bwa.
        TODO: Either manually add checking all binaries are downloaded properly
        or come up with automated way of testing this.
        """
        test_tools_dir = os.path.join(settings.TEMP_FILE_ROOT, 'test_tools')

        # # Samtools
        # download_tools(tools=['samtools'], tools_dir=test_tools_dir)
        # for binary_suffix in external_tool_registry.SAMTOOLS['bin_suffixes']:
        #     full_path_binary = binary_suffix % test_tools_dir
        #     assert os.path.exists(full_path_binary)

        # BWA
        download_tools(tools=['bwa'], tools_dir=test_tools_dir)
        for binary_suffix in external_tool_registry.BWA['bin_suffixes']:
            full_path_binary = binary_suffix % test_tools_dir
            assert os.path.exists(full_path_binary)
