import os
import pickle
import shutil
import sys
import tempfile

# Number of characters that snpeff requires per line of the genbank file
SNPEFF_MIN_LINE_LEN = 20

def setup_django_env():
    """Setting up the django/python env is a PITA but must be done."""
    # First make sure settings are on the path.
    PWD = os.path.dirname(os.path.realpath(__file__ ))
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


def fn_runner(fn, args_list, concurrent=False):
    """Helper method that handles calling a method depending on whether
    concurrent is True or not.

    Returns:
        If concurrent=True, immediately returns a celery.Result object without
        blocking. Otherwise blocks while executing the function, returning an
        implicit None.
    """
    from main.tasks import generic_task

    if concurrent:
        return generic_task.delay(fn.__name__, args_list)
    else:
        fn(*args_list)

import urllib2

def internet_on():
    try:
        response=urllib2.urlopen('http://74.125.228.100',timeout=1)
        return True
    except urllib2.URLError as err: pass
    return False

def ensure_line_lengths(gb_file):
    """
    Since SNPEFF requires at least 20 characters per line, we will add trailing
    spaces to every line of the genbank file and overwrite the file.
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
    

