import boto
import os
from boto.s3.connection import Key, S3Connection
from django.conf import settings
from models import Project
from functools import wraps
from contextlib import contextmanager
import tempfile
import hashlib
import logging


boto.set_stream_logger('boto')
logger = logging.getLogger('s3')
logging.getLogger('boto').setLevel(logging.INFO)

S3 = S3Connection(settings.AWS_SERVER_PUBLIC_KEY, settings.AWS_SERVER_SECRET_KEY)
aws_bucket = S3.get_bucket(settings.AWS_EXPECTED_BUCKET)


def calc_file_md5(filepath):
    md5 = hashlib.md5()
    with open(filepath,'rb') as f: 
        for chunk in iter(lambda: f.read(128*md5.block_size), b''): 
            md5.update(chunk)
    return md5.hexdigest()

def s3_get_directory(s3_dir, local_dir):
    rs = aws_bucket.list(s3_dir)
    logger.info("Getting s3://%s/%s to file://%s" % (
        aws_bucket.name, s3_dir, os.path.abspath(local_dir)))

    for key in rs:
        relpath = os.path.relpath(key.name, s3_dir)
        filepath = os.path.join(local_dir, relpath)
        directory = os.path.dirname(filepath)
        if not os.path.exists(directory):
            os.makedirs(directory)

        if os.path.isfile(filepath):
            if key.etag.strip("\"") == calc_file_md5(filepath):
                continue  # already fetched.
        key.get_contents_to_filename(filepath)

def s3_put_directory(s3_dir, local_dir):
    logger.info("Putting file://%s to s3://%s/%s" % (
        os.path.abspath(local_dir), aws_bucket.name, s3_dir))

    for dirname, dirnames, filenames in os.walk(local_dir):
        for filename in filenames:
            if filename == ".DS_Store":
                continue
            filepath = os.path.join(dirname, filename)
            relpath = os.path.relpath(filepath, local_dir)
            key = os.path.join(s3_dir, relpath)
            aws_key = aws_bucket.get_key(key)
            if aws_key:
                # assume the content of file did not change if md5 hashes are consistent.
                if aws_key.etag.strip("\"") == calc_file_md5(filepath):
                    continue
            else:
                aws_key = Key(aws_bucket, key)
            aws_key.set_contents_from_filename(filepath)
            aws_key.set_acl('public-read')

def s3_delete(key):
    aws_key = aws_bucket.get_key(key)
    aws_key.delete()

def s3_get_string(key):
    aws_key = aws_bucket.get_key(key)
    return aws_key.get_contents_as_string()

def s3_get(key, location):
    aws_key = aws_bucket.get_key(key)
    return aws_key.get_contents_to_filename(location)


@contextmanager
def s3_temp_get(s3file):
    aws_key = aws_bucket.get_key(s3file.key)
    assert aws_key

    f = tempfile.NamedTemporaryFile(delete=True, suffix=s3file.name)
    aws_key.get_contents_to_filename(f.name)

    yield f.name

    f.close()

@contextmanager
def project_s3_persisted(project):
    project = project
    assert isinstance(project, Project)

    if project.is_s3_backed():
        s3_get_directory(project.get_s3_model_data_dir(), project.get_model_data_dir())
    
    yield project

    if project.is_s3_backed():
        s3_put_directory(project.get_s3_model_data_dir(), project.get_model_data_dir())

def project_files_needed(func):
    """A decorator function to wrap function to make sure the availability and
    persistance of project files for the period of function execution.

    NOTE: 1st argument of function must be models.Project instance.
    """
    @wraps(func)
    def wrapper(project, *args, **kwargs):
        assert isinstance(project, Project)
        with project_s3_persisted(project):
            return func(project, *args, **kwargs)
    return wrapper
