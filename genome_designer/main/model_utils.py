"""
Utility objects and functions for interacting with models.
"""

import hashlib
from uuid import uuid4

from django.conf import settings
from django.db import IntegrityError
from django.db import models
from django.db import transaction

import os
import re
import stat


# Size of unique UIDs.
UUID_SIZE = 8

# Size of hash id for filenames that store long alt.
LONG_ALT_HASH_SIZE = 8


###############################################################################
# Mixins
###############################################################################

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


class JSONDataSubfieldsMixin(object):
    """Mixin that allows looking up attributes that are not directly on the
    model object to be looked for in JSON data.
    """

    def __getattr__(self, name):
        """Override.

        Automatically called if an attribute is not found in the typical
        place.

        Our implementation checks the data or metadata dict, raising an
        AttributeError if not found.

        See:
           http://docs.python.org/2/reference/datamodel.html#object.__getattr__
        """
        if hasattr(self, 'data'):
            data_attribute = self.data
        elif hasattr(self, 'metadata'):
            data_attribute = self.metadata
        else:
            raise AttributeError

        try:
            return data_attribute[name]
        except:
            raise AttributeError


###############################################################################
# Unique Uid
###############################################################################

class SafeCreateModelManager(models.Manager):
    """Manager that retries creation due to IntegrityError.
    """

    # The maximum number of clash resolutions allowed.
    # The probability of a clash for a # particular model (26 + 10) ^ 8 which
    # should not be a problem for now. If it takes more than 10 attempts to
    # resolve, either the db is too full, or there is a bug.
    MAX_UID_CLASHES = 10

    def create(self, uid_fail_count=0, **kwargs):
        """Create, allowing for rare uid clash.
        """
        try:
            return self.get_query_set().create(**kwargs)
        except IntegrityError as e:
            transaction.commit_unless_managed()
            uid_fail_count += 1
            if uid_fail_count > self.MAX_UID_CLASHES:
                raise e
            return self.create(uid_fail_count=uid_fail_count, **kwargs)


def short_uuid():
    """Generates a short uuid, repeating if necessary to maintain
    uniqueness.

    This is passed as the value of the default argument to any model
    with a uid field.

    NOTE: I thought I could make this a @classmethod but I don't think it's
    possible to access the instance class at the scope where model fields
    are declared.
    """
    return str(uuid4())[:UUID_SIZE]


class UniqueUidModelMixin(models.Model):
    """Mixin that adds a field for a unique uid and catches an IntegrityError
    and retries in case the same uid is used.

    NOTE: This overrides the manager. Further changes that require changing
    the manager will require considering the SafeCreateModelManager.
    """
    class Meta:
        abstract = True

    uid = models.CharField(max_length=UUID_SIZE, unique=True,
            default=short_uuid)

    # Override the default manager.
    objects = SafeCreateModelManager()


###############################################################################
# Misc helpers
###############################################################################

def ensure_exists_0775_dir(destination):
    """Creates a directory with 0775 permissions, and gets the group of the
    parent directory rather than the effective group of the process.

    The 0775 permissions are interpreted as all permissions for owner and
    group, with read-only permissions for all others.

    Does nothing if it already existed.  Errors if creation fails.
    """
    try:
        os.makedirs(destination)
        # User and group have all permissions | get group id from directory.
        os.chmod(destination, stat.S_ISGID | 0775)
    except OSError as e:
        if e.errno != 17: # File exists
            raise
        # Otherwise, nothing to do.
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


def get_dataset_with_type(entity, type, compressed=False):
    """Returns the Dataset with the requested type, or None if doesn't exist.

    If there are both compressed and uncompressed versions, return the 
    uncompressed unless the compressed is asked for.
    """
    results = [r for r in entity.dataset_set.filter(type=type) 
            if r.is_compressed() == compressed]

    assert len(results) < 2, ("More than one Datasets of type %s for entity %s."
            % (type, str(entity)))

    if len(results) > 0:
        return results[0]
    return None


def get_normalized_alt_representation(alt_str):
    """Returns a normalized represenation of the alt.

    This provides support for long alt values which we don't want to store in
    the database. Instead we store a unique hash that is the name of the file
    on disk.

    For a short alt <= 10bp, this is function just return the alt unchanged.

    For a long alt, this function returns a string of the form
    "LONG:<size>|FILE:<hash>", where <hash> is computed directly from alt_str.
    """
    if len(alt_str) <= 10:
        return alt_str

    size = len(alt_str)
    hash_part = hashlib.md5(alt_str).hexdigest()[:LONG_ALT_HASH_SIZE]
    return 'LONG:{size},HASH:{hash_part}'.format(
            size=size, hash_part=hash_part)


def get_long_alt_path(ref_genome, hash_part):
    """Returns path to file where long alt is stored.
    """
    long_alts_dir = ref_genome.get_long_variant_alts_dir()

    # NOTE: This line helps transitioning existing installations of Millstone
    # but it may be wonky with respect to concurrency.
    # TODO: Remove.
    ensure_exists_0775_dir(long_alts_dir)

    return os.path.join(long_alts_dir, hash_part + '.txt')
