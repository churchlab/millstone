"""
Custom model fields.
"""

import json

from django.db import models

from psycopg2.extras import Json as psycopg2_Json


class PostgresJsonField(models.Field):
    """Django representation of a json field.

    The current implementation is hard-coded to work with the Postgresl 9.3
    json field.

    The most popular implementation of a library that provides json support
    for Django, django-jsonfield (https://pypi.python.org/pypi/jsonfield), is
    actually currently available in pypi was designed to give a json
    abstraction in the Python code, but store the underlying value as a
    TextField equivalent. We specifically need the Postgresql 9.3
    json field so we implement it our own way here.
    """

    # Require for to_python() to work.
    # See: https://docs.djangoproject.com/en/1.5/howto/custom-model-fields/#the-subfieldbase-metaclass
    __metaclass__ = models.SubfieldBase

    def db_type(self, connection):
        return 'json'

    def to_python(self, value):
        """Converts the raw value returned from the database to a Python object.

        Right now, the default is to return a dictionary since everywhere
        we use this field we expect an object. A more general, and correct,
        implementation would be to intelligently deserialize or return None,
        supporting types other than dictionaries.
        """
        if isinstance(value, dict):
            return value
        if not len(value):
            return {}
        return json.loads(value)

    def get_db_prep_value(self, value, connection, prepared=False):
        """Converts value to a form that can be stored in a Potsgres json
        column.
        """
        if isinstance(value, dict):
            return psycopg2_Json(value)
        # Else assume it's json already.
        return value

    def south_field_triple(self):
        """Required for south to work.
        """
        from south.modelsinspector import introspector
        name = '%s.%s' % (self.__class__.__module__, self.__class__.__name__)
        args, kwargs = introspector(self)
        return name, args, kwargs
