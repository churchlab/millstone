"""
Miscellaneous utility functions.
"""
import sys
import os
import re
import shutil
import tempfile
import urllib2
import string

from Bio import SeqIO
from django.conf import settings

# Number of characters that snpeff requires per line of the genbank file
SNPEFF_MIN_LINE_LEN = 20


def make_temp_file(label, extension):
    if not os.path.exists(settings.TEMP_FILE_ROOT):
        os.mkdir(settings.TEMP_FILE_ROOT)
    _, temp_file_path = tempfile.mkstemp(
            suffix=('_' + generate_safe_filename_prefix_from_label(label) +
                    extension),
            dir=settings.TEMP_FILE_ROOT)
    return temp_file_path


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


def remove_whitespace(a_string):
    """Returns string with whitespace removed and stripped.
    """
    return re.sub('\W', '_', a_string.strip())


def uppercase_underscore(a_string):
    """
    Internally some strings that are Title Cased With Spaces should be
    UPPER_CASED_WITH_UNDERSCORES. 
    """
    a_string = a_string.replace(' ','_')
    return a_string.upper()


def lowercase_underscore(a_string):
    """Returns lowercased version of string,where whitespace is replaced by
    underscores and extra whitespace is stripped.
    """
    return re.sub('\W', '_', a_string.strip().lower())


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


def convert_fasta_to_fastq(fa_path, fq_path):
    """Generates a fasta file from a fastq file
    """
    with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [40] * len(record)
            SeqIO.write(record, fastq, "fastq")


def convert_seqrecord_to_fastq(seqrecord, fq_path):
    """Generates a fasta file from a Biopython seqrecord object
    """
    with open(fq_path, "w") as fastq:
        seqrecord.letter_annotations["phred_quality"] = [40] * len(seqrecord)
        SeqIO.write(seqrecord, fastq, "fastq")
