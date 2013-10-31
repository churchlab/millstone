import boto
from boto.s3.connection import Key, S3Connection
from django.conf import settings

boto.set_stream_logger('boto')
S3 = S3Connection(settings.AWS_SERVER_PUBLIC_KEY, settings.AWS_SERVER_SECRET_KEY)


def s3_delete(bucket, key):
    aws_bucket = S3.get_bucket(bucket, validate=False)
    aws_key = aws_bucket.get_key(key)
    aws_key.delete()

def s3_get(bucket, key, location):
    aws_bucket = S3.get_bucket(bucket)
    aws_key = aws_bucket.get_key(key)
    aws_key.get_contents_to_filename(location)

def s3_upload(bucket, key, location):
    aws_bucket = S3.get_bucket(bucket)
    aws_key = aws_bucket.get_key(key)
    aws_key.set_contents_from_filename(location)
