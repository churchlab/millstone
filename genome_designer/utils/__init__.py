"""
Miscellaneous utility functions.
"""

# __all__ = ['generate_safe_filename_prefix_from_label']

import os
import re
import shutil
import sys
import tempfile
import urllib2
import string

# Number of characters that snpeff requires per line of the genbank file
SNPEFF_MIN_LINE_LEN = 20


def setup_django_env():
    """Setting up the django/python env is a PITA but must be done."""
    # First make sure settings are on the path.
    PWD = os.path.dirname(os.path.realpath(__file__))
    settings_path = os.path.split(PWD)[0]
    os.sys.path.append(settings_path)

    # Now setup the environment using the settings
    from django.core import management
    try:
        import settings # Assumed to be in the os path
    except ImportError:
        sys.stderr.write(
                "Error: Can't find the file 'settings.py' in the directory " +
                "containing %r. It appears you've customized things.\nYou'll" +
                "have to run django-admin.py, passing it your settings " +
                "module.\n(If the file settings.py does indeed exist, it's " +
                "causing an ImportError somehow.)\n" % __file__)
        sys.exit(1)

    management.setup_environ(settings)

    # All of our scripts should maintain group writer permissions since we have
    # multiple different users potentially writing to files, though all are
    # in the genome-designer group.
    # We set the process umask here to enforce this.
    # TODO: Is there a better place to set this?
    os.umask(002)


def internet_on():
    """Check whether we're connected to the Internet.
    """
    try:
        urllib2.urlopen('http://74.125.228.100', timeout=1)
        return True
    except urllib2.URLError:
        pass
    return False


def ensure_line_lengths(gb_file):
    """Since SNPEFF requires at least 20 characters per line, we will add
    trailing spaces to every line of the genbank file and overwrite the file.
    """

    fh_tmp = tempfile.NamedTemporaryFile(delete=False)
    with open(gb_file, 'r') as fh_in:
        for line in fh_in:
            line = line.rstrip()
            if len(line) < SNPEFF_MIN_LINE_LEN:
                line += ' '*(SNPEFF_MIN_LINE_LEN - len(line))
            print >> fh_tmp, line
    fh_tmp.close()
    shutil.move(fh_tmp.name, gb_file)

def merge_nested_dictionaries(d1, d2, allow_update=True):
    '''
    Primarily for adding dictionaries to an existing json dict. 

    http://stackoverflow.com/
            questions/17857926/
            how-merge-with-appending-two-nested-dictionaries-in-python

    If update is true, keys will overwrite if not the same type.
    '''

    def merge_values(a,b):
        if a==None or b==None:
            return a or b
        # now handle cases where both have values
        if type(a) == dict and type(b) == dict:
            return add_dict(a, b)
        if type(a)==list and type(b) == list:
            return a + b
        if allow_update:
            return b
        else:
            raise Exception(
                    "Dicts don't match:\nA:{:s}\nB:{:s}".format(
                            a,b))

    def add_dict(d1,d2):
        return dict(
            [
                (key,
                 merge_values(
                     d1.get(key,None),
                     d2.get(key,None)))
                for key
                in set(d1.keys()).union(d2.keys())
            ])

    return merge_values(d1, d2)

def uppercase_underscore(a_string):
    """
    Internally some strings that are Title Cased With Spaces should be
    UPPER_CASED_WITH_UNDERSCORES. 
    """
    a_string = a_string.replace(' ','_')
    return a_string.upper()

def titlecase_spaces(a_string):
    """
    Externally some strings that are UPPER_CASED_WITH_UNDERSCORES should be
    Title Cased With Spaces.
    """
    a_string = a_string.replace('_',' ')
    return string.capwords(a_string)


def generate_safe_filename_prefix_from_label(label):
    """Generates safe filename prefix by avoiding non-alphanumeric chars.
    """
    return re.sub('\W', '_', label.lower())
