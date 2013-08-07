#!/usr/bin/env python

"""
Setup script.
"""

import os
import platform
import re
import shutil
import stat
import urllib
import zipfile

from bs4 import BeautifulSoup

from settings import JBROwSE_DATA_SYMLINK_PATH
from settings import MEDIA_ROOT
from settings import PWD as GD_ROOT
from settings import TOOLS_DIR


TOOLS_URLS = { 
    'bwa' : {
        'Darwin': ['https://www.dropbox.com/s/w6gs82gn0rmc60w/bwa'],
        'Linux': ['https://www.dropbox.com/s/eq3533wx8ay41sx/bwa'],
    },
    'freebayes' : {
        'Darwin' : [
            'https://www.dropbox.com/s/18us45ggflvci6e/bamleftalign',
            'https://www.dropbox.com/s/7pzjie4u4itxyjv/freebayes'
        ],
        'Linux' : [
            'https://www.dropbox.com/s/3qszohowh1gj1u0/bamleftalign',
            'https://www.dropbox.com/s/n2we6c6bmlbyi6f/freebayes'
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
            'https://www.dropbox.com/s/pjw0h0vzo2ls79t/samtools-0.1.19-linux.zip'
        ],
        'Darwin' : [
            'https://www.dropbox.com/s/vjtegtyncmq8a8o/samtools-0.1.19-darwin.zip'
        ]
    },
    'picard' : [
        'https://www.dropbox.com/s/189loxsmsch1xbb/picard-tools-1.96.zip'
    ]
}    

def setup():
    download_tools()
    setup_jbrowse()

def download_tools():
    # Create tools dir if it doesn't exist.
    if not os.path.isdir(TOOLS_DIR):
        os.mkdir(TOOLS_DIR)

    sys_type = platform.system()
    
    # Download the required tools from our Dropbox.
    for tool in TOOLS_URLS:
        
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
                tool_zipped = zipfile.ZipFile(tmp_path)
                tool_zipped.extractall(dest_dir)
                                
            # Otherwise download files to the correct location.
            else:
                urllib.urlretrieve(tool_url,dest_path)
                print '    %s' % dest_path
        
        # This is really hackish, but because ZipFile does not save file
        # permissions, we need to make the bins executable. Luckily,
        # in all cases (so far), the bin has the same name as the tool...
        tool_bin_file = os.path.join(dest_dir,tool)
        if os.path.exists(tool_bin_file):
            os.chmod(os.path.join(dest_dir,tool), stat.S_IXUSR)

def _get_file_url_from_dropbox(dropbox_url, filename):
    """
    Dropbox takes you to a download page when you share a link, so to get 
    past that, we use BeautifulSoup to get the file, as per goo.gl/EeOGjC
    """
    try:
        url_text = urllib.urlopen(dropbox_url).read()
        url_soup = BeautifulSoup(url_text)
    except:
        raise InputError(
            'Could not get file %s from Dropbox URL %s. Check settings.' % (
                filename, dropbox_url))

    return url_soup.find_all(href=re.compile(filename))[0]['href']


def setup_jbrowse():
    """Sets up JBrowse.

    TODO: Move remaining manual steps here.
    """
    if not os.path.exists(JBROwSE_DATA_SYMLINK_PATH):
        os.symlink(os.path.join(GD_ROOT, MEDIA_ROOT),
                JBROwSE_DATA_SYMLINK_PATH)


if __name__ == '__main__':
    setup()
