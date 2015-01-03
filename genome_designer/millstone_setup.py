#!/usr/bin/env python

"""
Setup script.
"""

import os
import platform
import re
import stat
import sys
import urllib
import zipfile

from settings import JBROWSE_DATA_SYMLINK_PATH
from settings import MEDIA_ROOT
from settings import PWD as GD_ROOT
from settings import TOOLS_DIR


TOOLS_URLS = {
    'bwa' : {
        'Darwin': ['https://www.dropbox.com/s/w6gs82gn0rmc60w/bwa'],
        'Linux': ['https://www.dropbox.com/s/eq3533wx8ay41sx/bwa'],
    },
    'freebayes' : {
        # we use the following tools in the freebayes git repo (+submodules):
        # * freebayes
        # * bamleftalign
        # * bamtools
        # * vcfstreamsort (from vcflib)
        # * vcfuniq (from vcflib)
        'Darwin' : [
            'https://www.dropbox.com/s/gjpqfr8n2601mfh/bamleftalign',
            'https://www.dropbox.com/s/jx2zdan4ci0mx6w/freebayes',
            'https://www.dropbox.com/s/f7twxy2r0h5ox20/bamtools',
            'https://www.dropbox.com/s/kvhkji81l16ukm5/vcfstreamsort',
            'https://www.dropbox.com/s/kaozc20jovxuzmk/vcfuniq'

        ],
        'Linux' : [
            'https://www.dropbox.com/s/3qszohowh1gj1u0/bamleftalign',
            'https://www.dropbox.com/s/n2we6c6bmlbyi6f/freebayes',
            'https://www.dropbox.com/s/5qrdi5prcjbz6js/bamtools',
            'https://www.dropbox.com/s/o0lpdbns3l94wy5/vcfstreamsort',
            'https://www.dropbox.com/s/j9y8ul4k71f4kx5/vcfuniq'
        ],
    },
    'snpEff' : [
        'https://www.dropbox.com/s/jnob4hkcwtquyag/snpEff-3.3e.zip',
        ],
    'gatk' : [
        #New GATK doesn't work for OSX 10.6 and below b/c it requires Java 7
        #'https://www.dropbox.com/s/47b0ivwlr80vtw1/GenomeAnalysisTK-2.6-5-gba531bd.zip'
        'https://www.dropbox.com/s/227curebr2prswn/GenomeAnalysisTK-legacy.zip'
        ],
    'samtools' : {
        'Linux' : [
            'https://www.dropbox.com/s/n3q9lpjgx224eix/samtools-0.1.20-linux.zip'
        ],
        'Darwin' : [
            'https://www.dropbox.com/s/wvd9iumzq1qi24h/samtools-0.1.20-darwin.zip'
        ]
    },
    'tabix' : {
        'Darwin' : [
            'https://www.dropbox.com/s/3mr5x3di2p513ma/tabix-0.2.4-darwin.zip'
        ],
        'Linux': [
            'https://www.dropbox.com/s/38y3p2e7za4conk/tabix-0.2.4-linux.zip'
        ]
    },
    'pindel' : {
        'Darwin' : [
            'https://www.dropbox.com/s/8wuh47i2hao3vlr/pindel-darwin.zip'
        ],
        'Linux' : [
            'https://www.dropbox.com/s/wyfve2dfgqbt94d/pindel-linux.zip'
        ]
    },
    'delly' : {
        'Darwin' : [
            'https://www.dropbox.com/s/l3qs0mbm2vapqf8/delly-darwin.zip'
        ],
        'Linux' : [
            'https://www.dropbox.com/s/k1gaxm1wb0odubc/delly-linux.zip'
        ]
    },
    'vcf-concat' : {
        'Darwin' : [
            'https://www.dropbox.com/s/ehowcdaic4tln0n/vcf-concat-darwin.zip'
        ],
        'Linux' : [
            'https://www.dropbox.com/s/buf80h0qfwjqaqp/vcf-concat-linux.zip'
        ]
    },
    'picard' : [
        'https://www.dropbox.com/s/189loxsmsch1xbb/picard-tools-1.96.zip'
    ],
    'lumpy' : {
        'Darwin' : [
            'https://www.dropbox.com/s/xyzmhtyjdfzhcza/lumpy-0.2.1-darwin.zip'
        ],
        'Linux' : [
            'https://www.dropbox.com/s/xl8honqq7gnh5yi/lumpy-0.2.1-linux.zip'
        ]
    },
    'fastqc' : [
            'https://www.dropbox.com/s/6d46rqjxyqi9k31/fastqc_v0.11.2.zip'
    ]
}


# For any tools that have multiple executables that need permission changes,
# for tools whose executables arent named after their tool
TOOLS_TO_EXECUTABLES = {
    'lumpy': ['*'], # make all executable
    'pindel': ['pindel', 'pindel2vcf'],
    'tabix': ['tabix', 'bgzip'],
    'freebayes': ['freebayes', 'bamtools', 'vcfstreamsort', 'vcfuniq']
}


def setup(arglist):
    if len(arglist) == 0:
        download_tools()
        setup_jbrowse()
    else:
        if 'jbrowse' in arglist:
            setup_jbrowse()
            arglist.remove('jbrowse')
        download_tools(tools=arglist)


def download_tools(tools=TOOLS_URLS.keys()):
    # Create tools dir if it doesn't exist.
    if not os.path.isdir(TOOLS_DIR):
        os.mkdir(TOOLS_DIR)

    sys_type = platform.system()

    # Download the required tools from our Dropbox.
    for tool in tools:

        print 'Grabbing %s from Dropbox...' % tool

        # If the tool is a dict, then it references different OSes (i.e.
        # compiled binaries). Download the correct one for your OS.
        if isinstance(TOOLS_URLS[tool],dict):
            try:
                tool_urls = TOOLS_URLS[tool][sys_type]
            except KeyError, e:
                raise(OSError(' '.join([
                    'Tool "%s" missing for your OS (%s),'
                    'OSs available are: (%s)']) % (
                        tool, sys_type, str(TOOLS_URLS[tool].keys()))))
        else:
            tool_urls = TOOLS_URLS[tool]

        # Make the directory for this tool, if it doesn't exist
        dest_dir = os.path.join(TOOLS_DIR,tool)

        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)

        for tool_url in tool_urls:

            # Grab last part of url after final '/' as file name
            dest_filename = re.match(r'^.*\/([^/]+)$', tool_url).group(1)
            dest_path = os.path.join(dest_dir,dest_filename)

            # Get the *actual* file URL from the dropbox popup
            tool_url = _get_file_url_from_dropbox(tool_url, dest_filename)

            # If *.zip, download to a tmp file and unzip it to the right dir.
            if dest_filename.endswith('.zip'):
                tmp_path, http_msg = urllib.urlretrieve(tool_url)
                print '    %s (Unzipping...)' % tmp_path
                try:
                    tool_zipped = zipfile.ZipFile(tmp_path)
                    tool_zipped.extractall(dest_dir)
                except:
                    'File {} missing or invalid format!'.format(
                            tool_url)

            # Otherwise download files to the correct location.
            else:
                urllib.urlretrieve(tool_url,dest_path)
                print '    %s' % dest_path

        # This is really hackish, but because ZipFile does not save file
        # permissions, we need to make the bins executable. In cases where
        # the executable name is different than the tool key in TOOL_URLS,
        # use TOOLS_TO_EXECUTABLES.
        if not tool in TOOLS_TO_EXECUTABLES:
            tool_bin_files = [tool]
        else:
            tool_bin_files = [exe for exe in TOOLS_TO_EXECUTABLES[tool]]
            if '*' in tool_bin_files:
                # All files are executables.
                tool_bin_files = [exe for exe in os.listdir(dest_dir)]

        for tbf in tool_bin_files:
            tool_bin_path = os.path.join(dest_dir, tbf)
            if os.path.exists(tool_bin_path):
                os.chmod(tool_bin_path,
                        stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)


def _get_file_url_from_dropbox(dropbox_url, filename):
    """Dropbox now supports modifying the shareable url with a simple
    param that will allow the tool to start downloading immediately.
    """
    return dropbox_url + '?dl=1'


def setup_jbrowse():
    """Sets up JBrowse.
    """
    # Unlink if there was a link and then create a new link.
    try:
        os.unlink(JBROWSE_DATA_SYMLINK_PATH)
    except OSError:
        # There was no symlink. That's fine.
        pass

    # Re-create the link.
    os.symlink(os.path.join(GD_ROOT, MEDIA_ROOT),
            JBROWSE_DATA_SYMLINK_PATH)

if __name__ == '__main__':
    setup(sys.argv[1:])
