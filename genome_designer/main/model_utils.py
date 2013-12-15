"""
Utility objects and functions for interacting with models.
"""

from uuid import uuid4

from django.conf import settings

import os
import pickle
import re
import stat

class VisibleFieldMixin(object):
    """Mixin that provides methods for controlling which fields are displayed.

    NOTE(gleb): I'm not sure if this is general enough to apply to other
        models.
    """
    @classmethod
    def get_field_order(clazz, **kwargs):
        if not hasattr(clazz, 'default_view_fields'):
            return []

        if not 'additional_field_list' in kwargs:
            return clazz.default_view_fields()

        # Otherwise, incorporate additional fields, without repetition.
        default_field_names = set([field_obj['field'] for field_obj in
                clazz.default_view_fields()])
        additional_fields = [
                {'field': field_name} for field_name
                in kwargs['additional_field_list']
                if not field_name in default_field_names]
        return clazz.default_view_fields() + additional_fields


def short_uuid(cls):
    """Generates a short uuid, repeating if necessary to maintain
    uniqueness.

    This is passed as the value of the default argument to any model
    with a uid field.

    NOTE: I thought I could make this a @classmethod but I don't think it's
    possible to access the instance class at the scope where model fields
    are declared.
    """
    UUID_SIZE = 8

    # Even with a short id, the probability of collision is very low,
    # but timeout just in case rather than risk locking up.
    timeout = 0
    while timeout < 1000:
        initial_long = str(uuid4())
        candidate = initial_long[:UUID_SIZE]
        if len(cls.objects.filter(uid=candidate)) == 0:
            return candidate
        else:
            timeout += 1
    raise RuntimeError, "Too many short_uuid attempts."


def ensure_exists_0775_dir(destination):
    """Creates a directory with 0775 permissions, and gets the group of the
    parent directory rather than the effective group of the process.

    The 0775 permissions are interpreted as all permissions for owner and
    group, with read-only permissions for all others.

    Does nothing if it already existed.  Errors if creation fails.
    """
    if not os.path.exists(destination):
        os.makedirs(destination)
        # User and group have all permissions | get group id from directory.
        os.chmod(destination, stat.S_ISGID | 0775)

    return True


def make_choices_tuple(type_class):
    """Creates a tuple of tuples object used as the choices attribute for
    a CharField that we want to limit choices.

    Args:
        type_class: A class with the attributes defining the types.
    """
    return tuple([
            (type_name, type_name) for type_name in dir(type_class)
            if not re.match(r'__*', type_name)
    ])


def assert_unique_types(type_class):
    """Function called at runtime to make sure types are unique.
    """
    all_type_name_list = [type_name for type_name in dir(type_class)
            if not re.match(r'__*', type_name)]
    assert len(all_type_name_list) == len(set(all_type_name_list))


def clean_filesystem_location(filesystem_location):
    """If the filesystem location contains the full absolute path,
    trim it to be relative to the Django app MEDIA_ROOT.
    """
    clean_filesystem_location = filesystem_location
    match = re.search(settings.MEDIA_ROOT, filesystem_location)
    if match:
        clean_filesystem_location = clean_filesystem_location[match.end() + 1:]
    return clean_filesystem_location


def get_dataset_with_type(entity, type):
    """Returns the Dataset with the requested type, or None if doesn't exist.
    """
    results = entity.dataset_set.filter(type=type)
    assert len(results) < 2, ("More than 2 Datasets of type %s for entity %s."
            % (str(entity), type))
    if len(results) > 0:
        return results[0]
    return None


def auto_generate_short_name(long_name):
    """Helper method to compute a short name from a long name."""
    SHORT_NAME_CHARS = 12
    tokens = [token.lower() for token in long_name.split()]
    short_name = '_'.join(tokens)
    short_name = short_name[:SHORT_NAME_CHARS]
    return short_name


def get_flattened_unpickled_data(data):
    """Returns a dictionary from key to string values.

    Tries to unpickle the values if possible.
    """
    clean_data = {}
    for key, value in data.iteritems():
        try:
            clean_value = pickle.loads(str(value))
        except:
            clean_value = str(value)
        clean_data[key] = clean_value
    return clean_data
