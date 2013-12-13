"""Database models.

Note on filesystem directory structure: (IN PROGRESS)

    Since we are storing data output from various bioinformatics programs, the
    models below result in the creation and maintenance of a directory
    structure for data location. In general, the strategy is to have a
    directory corresponding to a model instance when possible. Hierarchy
    is used when models can be properly nested.

    An example layout for a user's data might look like:
        ../projects/1324abcd/
        ../projects/1324abcd/alignments/
        ../projects/1324abcd/samples/
        ../projects/1324abcd/samples/1234abcd
        ../projects/1324abcd/samples/5678jklm
        ../projects/1324abcd/ref_genomes/
        ../projects/1324abcd/variant_calls/

Implementation Notes:
    * get_field_order() for each model/table is used by the Adapter class
      in adapters.py to know WHICH FIELDS are to be displayed and WHAT ORDER.
      If you don't return a field in get_field_order, it won't be sent to
      datatables.js for display.
      Each field consists of a dict with a 'field' key, which is the name of
      the field, and an optional 'verbose' key, which is the display name of
      the field in the datatable. If 'verbose' is absent, then the underscores
      are converted to spaces and each word is Title Cased.

"""

import os
import pickle
import re
import shutil
import stat
from uuid import uuid4
from contextlib import contextmanager
import tempfile

from django.conf import settings
from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.db import models
from django.db.models import Model
from django.db.models.signals import post_delete
from django.db.models.signals import post_save
from jsonfield import JSONField


###############################################################################
# Utilities
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


###############################################################################
# User-related models
###############################################################################

class UserProfile(Model):
    """A UserProfile which is separate from the django auth User.

    This references the auth.User and opens up the possibility of
    adding additional fields.
    """
    # A one-to-one mapping to the django User model.
    user = models.OneToOneField(User)

    # A unique id, for urls or filesystem locations.
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(UserProfile)))

    def __unicode__(self):
        return self.user.username

# Since the registration flow creates a django User object, we want to make
# sure that the corresponding UserProfile is also created
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        user_profile = UserProfile.objects.create(user=instance)

post_save.connect(create_user_profile, sender=User,
        dispatch_uid='user_profile_create')


###############################################################################
# Data wrappers
###############################################################################

class Dataset(Model):
    """A specific data file with a location on the filesystem.

    Basically a wrapper for a file on the file system.

    This is similar to the Galaxy notion of a dataset.
    """
    # NOTE: I'm not sure whether we'll need uid for this, but keeping it
    # just in case for now.
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(Dataset)))

    # The type of data this represents (e.g. Dataset.Type.BWA_ALIGN).
    # This is a semantic identifier for the kinds of operations
    # that can be performed with this Dataset.
    class TYPE:
        """The type of this dataset.

        Limit to 40-chars as per Dataset.type field def.
        """
        REFERENCE_GENOME_FASTA = 'rgf' # fasta
        REFERENCE_GENOME_GENBANK = 'rgg' #genbank
        FASTQ1 = 'f1'
        FASTQ2 = 'f2'
        BWA_ALIGN = 'bwa_align'
        BWA_ALIGN_ERROR = 'bwa_align_error'
        VCF_FREEBAYES = 'vcff'
        VCF_USERINPUT = 'vcfu'
        VCF_FREEBAYES_SNPEFF = 'vcffe'

    TYPE_CHOICES = make_choices_tuple(TYPE)
    type = models.CharField(max_length=40, choices=TYPE_CHOICES)

    # Human-readable identifier. Also used for JBrowse.
    label = models.CharField(max_length=256)

    # Location on the filesystem relative to settings.MEDIA_ROOT.
    filesystem_location = models.CharField(max_length=512, blank=True)

    # When the dataset is a result of a computation, we'll set a status on it.
    # NOTE: The reliability of the present implementation of this model feature
    # is questionable.
    class STATUS:
        """The status of running this Dataset.

        Limit to 40-chars as per Dataset.status field def.
        """
        UNKNOWN = 'UNKNOWN'
        COMPUTING = 'COMPUTING'
        READY = 'READY'
        FAILED = 'FAILED'
    STATUS_CHOICES = make_choices_tuple(STATUS)
    status = models.CharField(max_length=40, choices=STATUS_CHOICES,
            default=STATUS.READY)

    # Dictionary of compression suffixes and programs to use to decompress
    # to a pipe
    COMPRESSION_TYPES = {
        '.gz': ('gzip', '-dc'),
        '.bz2': ('bzcat',),
        '.zip': ('unzip','-p')
    }

    def __unicode__(self):
        return self.label

    def get_absolute_location(self):
        """Returns the full path to the file on the filesystem.
        """
        return os.path.join(settings.PWD, settings.MEDIA_ROOT,
                self.filesystem_location)

    def is_compressed(self):
        """
        Checks file path for .bz2 or .gz ending, and if so, returns true.
        """
        return self.filesystem_location.endswith(
                tuple(self.COMPRESSION_TYPES.keys()))

    def wrap_if_compressed(self):
        """ This helper function returns a process substitution string
        to be used by /bin/bash if the fastq read file is compressed, otherwise
        it just returns get_absolute_location().
        """
        absolute_location = self.get_absolute_location()
        if self.is_compressed():
            extension = os.path.splitext(self.filesystem_location)[1]
            program = ' '.join(self.COMPRESSION_TYPES[extension])
            return '<({:s} {:s})'.format(
                    program, absolute_location)
        else:
            return absolute_location

    @contextmanager
    def stream(self):
        """
        If dataset is compressed, return a named pipe that decompressed the
        file, else just return the absolute location.
        """

        raise NotImplementedError
        # Currently the below isn't working; the mkfifo blocks itself so I can't
        # seem to read and write to it at the same time. For now, we've decided
        # to go for process substitution in Bash (see wrap_if_compressed(),
        # although this requires the use of Shell=True.

        # dirname = tempfile.mkdtemp()
        # p = None

        # try:
        #     if not self.is_compressed():
        #         return self.get_absolute_location()

        #     path = os.path.join(dirname, 'named_pipe')
        #     os.mkfifo(path)
        #     extension = os.path.splitext(self.filesystem_location)[1]
        #     program = self.COMPRESSION_TYPES[extension]
        #     with open(path, 'w') as wpipe:
        #         p = Popen(program.append(path)) # write to path
        #         return path
        # finally:
        #     shutil.rmtree(dirname)
        #     if p: p.close()

# Make sure the Dataset types are unique. This runs once at startup.
assert_unique_types(Dataset.TYPE)

###############################################################################
# Project models
###############################################################################

class Project(Model):
    """A single project belonging to a user.

    A project groups together ReferenceGenomes, ExperimentSamples, and other
    data generated by tools during analysis.
    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(Project)))

    # The project owner.
    # TODO: Implement permissions system so that projects can be shared.
    owner = models.ForeignKey('UserProfile')

    # The human-readable title of the project.
    title = models.CharField(max_length=256)

    s3_backed = models.BooleanField(default=settings.S3_ENABLED)

    def __unicode__(self):
        return self.title + '-' + str(self.owner)

    def is_s3_backed(self):
        return self.s3_backed

    def get_s3_model_data_dir(self):
        return os.path.join("projects", str(self.uid))

    def get_model_data_root(self):
        """Get the absolute location where all project data is stored.
        """
        return os.path.join(settings.PWD, settings.MEDIA_ROOT, 'projects')

    def get_model_data_dir(self):
        """Get the full path to where the user's data is stored.

        The data dir is the media root url combined with the user id.
        """
        return os.path.join(self.get_model_data_root(), str(self.uid))

    def ensure_model_data_dir_exists(self):
        """Ensure that a data directory exists for the user.

        The data directory is named according to the UserProfile.id.
        """
        # Make sure the root of projects exists
        ensure_exists_0775_dir(self.get_model_data_root())

        # Check whether the data dir exists, and create it if not.
        return ensure_exists_0775_dir(self.get_model_data_dir())

    def delete_model_data_dir(self):
        """Removes all data associated with this model.

        WARNING: Be careful with this method!
        """
        data_dir = self.get_model_data_dir()
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'uid'},
                {'field':'title'}]


# When a new Project is created, create the data directory.
def post_project_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_project_create, sender=Project,
        dispatch_uid='project_create')

# Delete all Project data when it is deleted.
def post_project_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
post_delete.connect(post_project_delete, sender=Project,
        dispatch_uid='project_delete')


class ReferenceGenome(Model):
    """A reference genome relative to which alignments and analysis are
    performed.
    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(ReferenceGenome)))

    # A ReferenceGenome belongs to a single Project.
    project = models.ForeignKey('Project')

    # A human-readable label for this genome.
    label = models.CharField(verbose_name="Name", max_length=256)

    # Number of chromosomes.
    num_chromosomes = models.IntegerField(verbose_name="# Chromosomes")

    # Number of chromosomes.
    num_bases = models.BigIntegerField(verbose_name="Total Size")

    # Datasets pointing to files on the system (e.g. .fasta files, etc.)
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
        verbose_name="Datasets")

    # a key/value list of all possible VCF fields, stored as a JSONField
    # and dynamically updated by dynamic_snp_filter_key_map.py
    variant_key_map = JSONField()

    def __unicode__(self):
        return self.label

    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'genome_designer.main.views.reference_genome_view',
                args=(self.project.uid, self.uid))

    def get_model_data_root(self):
        """Get the root location for all data of this type in the project.
        """
        return os.path.join(self.project.get_model_data_dir(), 'ref_genomes')

    def get_model_data_dir(self):
        """Get the full path to the location of this model's data.
        """
        return os.path.join(self.get_model_data_root(), str(self.uid))

    def ensure_model_data_dir_exists(self):
        """Ensure that a data directory exists for this model.
        """
        # Make sure the root exists.
        ensure_exists_0775_dir(self.get_model_data_root())

        # Check whether the data dir exists, and create it if not.
        return ensure_exists_0775_dir(self.get_model_data_dir())

    def get_jbrowse_directory_path(self):
        """Returns the full path to the root of JBrowse data for this
        ReferenceGenome.
        """
        return os.path.join(self.get_model_data_dir(), 'jbrowse')

    def ensure_jbrowse_dir(self):
        """Ensures that the jbrowse data dir exists."""
        return ensure_exists_0775_dir(self.get_jbrowse_directory_path())

    def get_snpeff_directory_path(self):
        """Returns the full path to the root of snpeff data for this
        ReferenceGenome.
        """
        return os.path.join(self.get_model_data_dir(), 'snpeff',
                self.uid)

    def ensure_snpeff_dir(self):
        """Ensures that the snpeff data dir exists."""
        return ensure_exists_0775_dir(self.get_snpeff_directory_path())

    def get_client_jbrowse_data_path(self):
        if self.project.is_s3_backed():
            return os.path.join(
                    'http://%s.s3.amazonaws.com/' % settings.S3_BUCKET,
                    'projects',
                    str(self.project.uid),
                    'ref_genomes',
                    str(self.uid),
                    'jbrowse')
        else:
            # Allow forcing through nginx (dev only).
            maybe_force_nginx = ''
            if settings.DEBUG_FORCE_JBROWSE_NGINX:
                maybe_force_nginx = 'http://localhost'

            return os.path.join(
                    maybe_force_nginx + '/jbrowse/index.html?data=gd_data/',
                    'projects',
                    str(self.project.uid),
                    'ref_genomes',
                    str(self.uid),
                    'jbrowse')

    def get_client_jbrowse_link(self):
        """Returns the link to jbrowse for this ReferenceGenome.

        Example url for user with uid 'abc', and project id 'xyz', and
        refgenome id 456:
            '/jbrowse/?data=gd_data/abc/projects/xyz/ref_genomes/456/jbrowse/'
        """
        return '/jbrowse/index.html?data=' + self.get_client_jbrowse_data_path()

    def is_annotated(self):
        """For several steps (notably snpEff), we want to check that this
        ReferenceGenome is annotated (i.e. it has a genbank file associated
        with it). This function returns True if a genbank file is available.
        """
        return self.dataset_set.filter(
                type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).exists()

    def get_variant_caller_common_map(self):
        return self.variant_key_map['snp_caller_common_data']

    def get_variant_alternate_map(self):
        return self.variant_key_map['snp_alternate_data']

    def get_variant_evidence_map(self):
        return self.variant_key_map['snp_evidence_data']

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'label'},
                {'field':'num_chromosomes', 'verbose':'# Chromosomes'},
                {'field':'num_bases', 'verbose':'Total Size'}]

class ExperimentSample(Model):
    """Model representing data for a particular experiment sample.

    Usually this corresponds to a pair of fastq reads for a particular 
    bacterial clone or colony, after barcode removal/de-multiplexing.
    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(ExperimentSample)))

    # A Sample belongs to a single Project.
    project = models.ForeignKey('Project')

    # Human-readable identifier.
    label = models.CharField('Sample Name', max_length=256)

    # Human-readable sample group that this value is in.
    group = models.CharField('Plate/Group', max_length=256)

    # Human-readable 'position' (well number, etc) that this sample is in
    # within a group
    well = models.CharField('Position/Well', max_length=256)

    # Number of reads in the sample.
    num_reads = models.BigIntegerField('# Reads', default=-1)

    # The datasets associated with this sample. The semantic sense of the
    # dataset can be determined from the Dataset.type field.
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
        verbose_name="Datasets")

    def __unicode__(self):
        return self.label

    def get_model_data_root(self):
        """Get the root location for all data of this type in the project.
        """
        return os.path.join(self.project.get_model_data_dir(), 'samples')

    def get_model_data_dir(self):
        """Get the full path to the location of this model's data.
        """
        return os.path.join(self.get_model_data_root(), str(self.uid))

    def ensure_model_data_dir_exists(self):
        """Ensure that a data directory exists for this model.
        """
        # Make sure the root of projects exists
        ensure_exists_0775_dir(self.get_model_data_root())

        # Check whether the data dir exists, and create it if not.
        return ensure_exists_0775_dir(self.get_model_data_dir())

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'label'},
                {'field':'group'},
                {'field':'well'}]


class AlignmentGroup(Model):
    """Collection of alignments of several related ExperimentSamples to the
    same ReferenceGenome.

    The reason for grouping alignments together is that our variant operations
    are generally relative to a single reference genome, and further, variant
    calling tools often take multiple alignments as input, thus it makes sense
    to group these in the database.

    For a one-to-one mapping of Alignment to Sample, see
    ExperimentSampleToAlignment.
    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(AlignmentGroup)))

    # Human-readable identifier.
    label = models.CharField(max_length=256, blank=True)

    # All alignments in this set are relative to this genome.
    reference_genome = models.ForeignKey('ReferenceGenome')

    # The aligner tool used for this alignment.
    class ALIGNER:
        """Constants for representing the aligner type.
        """
        BWA = 'BWA'
    ALIGNER_CHOICES = make_choices_tuple(ALIGNER)
    aligner = models.CharField(max_length=10, choices=ALIGNER_CHOICES)

    # Times for the alignment run.
    start_time = models.DateTimeField(auto_now=True)
    end_time = models.DateTimeField(auto_now=True)

    # Datasets pointing to files on the system (e.g. .fasta files, etc.)
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
            verbose_name="Datasets")

    @property
    def status(self):
        """Returns an aggregate status.
        """
        any_computing = False
        for alignment in self.experimentsampletoalignment_set.all():
            alignment_dataset = get_dataset_with_type(alignment,
                    Dataset.TYPE.BWA_ALIGN)
            if not alignment_dataset:
                return Dataset.STATUS.UNKNOWN
            if (alignment_dataset.status in [
                    Dataset.STATUS.FAILED, Dataset.STATUS.UNKNOWN]):
                return alignment_dataset.status
            if alignment_dataset.status == Dataset.STATUS.COMPUTING:
                any_computing = True
        if any_computing:
            return Dataset.STATUS.COMPUTING
        return Dataset.STATUS.READY

    def __unicode__(self):
        return self.label

    def get_model_data_root(self):
        """Get the root location for all data of this type in the project.
        """
        return os.path.join(self.reference_genome.project.get_model_data_dir(),
                'alignment_groups')

    def get_model_data_dir(self):
        """Get the full path to the location of this model's data.
        """
        return os.path.join(self.get_model_data_root(), str(self.uid))

    def ensure_model_data_dir_exists(self):
        """Ensure that a data directory exists for this model.
        """
        # Make sure the root exists.
        ensure_exists_0775_dir(self.get_model_data_root())

        # Check whether the data dir exists, and create it if not.
        return ensure_exists_0775_dir(self.get_model_data_dir())

    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'genome_designer.main.views.alignment_view',
                args=(self.reference_genome.project.uid, self.uid))

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'label'},
                {'field':'reference_genome'},
                {'field':'aligner'},
                {'field':'status', 'verbose':'Job Status'},
                {'field':'start_time'},
                {'field':'end_time'}]


class ExperimentSampleToAlignment(Model):
    """Model that describes the alignment of a single ExperimentSample
    to an AlignmentGroup.
    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(ExperimentSampleToAlignment)))

    alignment_group = models.ForeignKey('AlignmentGroup')

    experiment_sample = models.ForeignKey('ExperimentSample')

    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True)

    @property
    def status(self):
        """The status of a running alignment job.
        """
        alignment_datasets = self.dataset_set.filter(
                type=Dataset.TYPE.BWA_ALIGN)
        if len(alignment_datasets) > 0:
            return alignment_datasets[0].status
        return 'UNDEFINED'

    @property
    def error_link(self):
        return ('<a href="' +
                reverse(
                        'genome_designer.main.views.sample_alignment_error_view',
                        args=(self.alignment_group.reference_genome.project.uid,
                                self.alignment_group.uid,
                                self.uid)) +
                        '">error</a>')

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [
            {'field':'experiment_sample'},
            {'field':'status', 'verbose':'Job Status'},
            {'field':'error_link', 'verbose': 'Error output', 'is_href': True},
        ]

    def get_model_data_root(self):
        """Get the root location for all data of this type in the project.
        """
        return os.path.join(self.alignment_group.get_model_data_dir(),
                'sample_alignments')

    def get_model_data_dir(self):
        """Get the full path to the location of this model's data.
        """
        return os.path.join(self.get_model_data_root(), str(self.uid))


    def ensure_model_data_dir_exists(self):
        """Ensure that a data directory exists for this model.
        """
        # Make sure the root exists.
        ensure_exists_0775_dir(self.get_model_data_root())

        # Check whether the data dir exists, and create it if not.
        return ensure_exists_0775_dir(self.get_model_data_dir())


    def delete_model_data_dir(self):
        """Removes all data associated with this model.

        WARNING: Be careful with this method!
        """
        data_dir = self.get_model_data_dir()
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)


# Delete all alignment data when the model is deleted.
def post_sample_alignment_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
post_delete.connect(post_sample_alignment_delete,
        sender=ExperimentSampleToAlignment,
        dispatch_uid='sample_alignment_delete')


###############################################################################
# Variants (SNVs and SVs)
###############################################################################

class Variant(Model):
    """An instance of a variation relative to a reference genome.

    This might be, for example, a SNV (single nucleotide variation) or a bigger
    SV (structural variation). We are intentionally using a unified model for
    these two classes of variations as the fundamental unit of genome analysis
    is really a diff.

    TODO: See code from Gemini paper (PLOS ONE 7/18/13) for ideas.

    A variant need not necessarily be associated with a specific sample; the
    VariantToExperimentSample model handles this association.
    """
    # Maybe used in url for view of this entity.
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(Variant)))

    class TYPE:
        DELETION = 'DELETION'
        INSERTION = 'INSERTION'
        TRANSITION = 'TRANSITION'
        TRANSVERSION = 'TRANSVERSION'
        COMPLEX = 'COMPLEX' # Multi-base in different genomes
    TYPE_CHOICES = make_choices_tuple(TYPE)
    type = models.CharField('Type', max_length=40, choices=TYPE_CHOICES)

    reference_genome = models.ForeignKey('ReferenceGenome',
        verbose_name='Reference Genome')

    chromosome = models.CharField('Chromosome', max_length=256, blank=True)

    position = models.BigIntegerField('Position')

    ref_value = models.TextField('Ref')

    def __init__(self, *args, **kwargs): 
        """If we are passed an alt_value field, we need to get_or_create
        VariantAlternate objects corresponding to them, and link them  up to
        this new variant. We're ignoring the handling the rare situation when a
        Variant has no alt_values, which we don't really want to happen. It is
        difficult to handle because sometimes the VariantAlternate objects are
        declared separately and added to the Variant after __init__()."""

        alts = kwargs.get('alt_value', None)

        # Here I'm mutating kwargs to get rid of alt_value, but I can't think
        # of a reason why this would be a problem, since we've already saved it.
        kwargs.pop('alt_value',None)

        # call super's __init__ without the alt_value field if present
        super(Variant, self).__init__(*args, **kwargs)

        if alts is None: return

        #handle case of one or multiple alt_values
        if not isinstance(alts, basestring):
            # alt_value is a list of alts
            alt_values = alts
        else:
            # alt value is one alt
            alt_values = [alts]

        for alt_value in alt_values:
            self.variantalternate_set.add(
                    VariantAlternate.objects.create(
                            variant=self,
                            alt_value=alt_value
                    )
            )

    @property
    def label(self):
        # NOTE: If adding a new VCCD object to a variant, this could change by
        # the addition of new variants. Is that an issue?
        return (
                str(self.position) + 
                '_' + self.ref_value + 
                '_' + ','.join(self.get_alternates()))

    def get_alternates(self):
        """ Return a base string for each alternate for this variant. """

        return [alt.alt_value for alt in self.variantalternate_set.all()]

    @property
    def jbrowse_link(self):
        ref_genome_jbrowse = self.reference_genome.get_client_jbrowse_link()
        location_param = '&loc=' + str(self.position)
        full_href = ref_genome_jbrowse + location_param
        return '<a href="' + full_href + '">jbrowse</a>'


    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'genome_designer.main.views.single_variant_view',
                args=(self.reference_genome.project.uid,
                        self.reference_genome.uid, self.uid))

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'label'},
                {'field':'jbrowse_link', 'verbose': 'JBrowse'},
                {'field':'chromosome'},
                {'field':'position'},
                # This is useless right now, always 'UNKNOWN'
                #{'field':'type'},
                {'field':'ref_value', 'verbose':'Ref'},
                {'field':'variantset_set',
                    'verbose':'Set Membership',
                    'classes':['label']}]



class VariantCallerCommonData(Model, VisibleFieldMixin):
    """Model that describes data provided by a specific caller about a
    particular Variant.

    The reason for this model is that most variant callers are run for multiple
    ExperientSamples at the same time, generating some common data for each
    variant found, as well as data unique to each ExperimentSample. This model
    represents the common shared data.

    To be even more specific, the VCF format typically gives a row for each
    variant, where some of the columns describe the variant, while other
    columns have a one-to-one correspondence to the alignments provided.
    """
    # Variant this object refers to. It's possible for multiple callers report
    # the same Variant so this is a many-to-one relationship.
    variant = models.ForeignKey('Variant')

    # Source dataset for this data.
    source_dataset = models.ForeignKey('Dataset')

    # Catch-all key-value data store.
    data = JSONField()

    def __getattr__(self, name):
        """Automatically called if an attribute is not found in the typical
        place.

        Our implementation checks the data dict, return the string 'undefined'
        if the value is not found.

        NOTE: Typically __getattr__ should raise an AttributeError if the value
        cannot be found, but the noisy nature or our data means returning
        'undefined' is more correct.

        See: http://docs.python.org/2/reference/datamodel.html#object.__getattr__
        """
        try:
            return pickle.loads(self.data[name])
        except:
            raise AttributeError

    def as_dict(self):
        """Converts a common data object into a dictionary from key to cleaned
        values.

        Cleaned generally means that fields that had to be pickled for storage
        are unpickled.

        Note that the original object had some SQL-level fields, but most of the
        data pased from the vcf file is in a catch-all 'data' field.
        This method flattens the structure so that all data is available on the
        resulting top-level object.

        Returns:
            A flattened dictionary of cleaned data.
        """
        return get_flattened_unpickled_data(self.data)

    @classmethod
    def default_view_fields(clazz):
        return []


class VariantAlternate(Model, VisibleFieldMixin):
    """A model listing alternate alleles for each variant."""

    uid = models.CharField(max_length=36,
        default=(lambda: short_uuid(VariantAlternate)))
    
    # Null is true here because we are adding this relationship during Variant's
    # overloaded __init__() so it hasn't been saved() yet. Otherwise it throws
    # an django.db.utils.IntegrityError: 
    #    main_variantalternate.variant_id may not be NULL
    variant = models.ForeignKey('Variant', null=True)

    alt_value = models.TextField('Alt')

    is_primary = models.BooleanField(default='False')

    # this json fields holds all PER ALT data (INFO data with num -1)
    data = JSONField()

    def __unicode__(self):
        return 'var: ' + str(self.variant) + ', alt:' + self.alt_value

    def as_dict(self):
        """Converts a alternate object into a dictionary from key to cleaned
        values.

        Cleaned generally means that fields that had to be pickled for storage
        are unpickled.

        Note that the original object had some SQL-level fields, but most of the
        data pased from the vcf file is in a catch-all 'data' field.
        This method flattens the structure so that all data is available on the
        resulting top-level object.

        Returns:
            A flattened dictionary of cleaned data.
        """
        return get_flattened_unpickled_data(self.data)

    # TODO: Do we want to explicitly link each VariantAlternate to
    # it's variant index in each VCCD object or VE object?
    # Currently it's done implicitly through the VCCD's data['ALT']
    # field and VE's data['gt_bases'] and data['GT'] fields, but these
    # are not checked for consistency. 
    @classmethod
    def default_view_fields(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'alt_value', 'verbose':'Alt(s)'}]

    
class VariantEvidence(Model, VisibleFieldMixin):
    """
    Evidence for a particular variant occurring in a particular
    ExperimentSample.
    """
    # Maybe used in url for view of this entity.
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(VariantEvidence)))

    # The specific ExperimentSample that this object provides evidence
    # of the respective variant occurring in.
    # NOTE: This implies the ReferenceGenome.
    experiment_sample = models.ForeignKey('ExperimentSample')

    # The location of the common data for this call.
    variant_caller_common_data = models.ForeignKey('VariantCallerCommonData')

    # One or more alternate alleles for this variant - 
    # Multiple are possible if the allele is called for multiple alts
    variantalternate_set = models.ManyToManyField('VariantAlternate')

    # Catch-all key-value set of data.
    # TODO: Extract interesting keys (e.g. gt_type) into their own SQL fields.
    data = JSONField()

    def __getattr__(self, name):
        """Automatically called if an attribute is not found in the typical
        place.

        Our implementation checks the data dict, return the string 'undefined'
        if the value is not found.

        NOTE: Typically __getattr__ should raise an AttributeError if the value
        cannot be found, but the noisy nature or our data means returning
        'undefined' is more correct.

        See: 
           http://docs.python.org/2/reference/datamodel.html#object.__getattr__

        """
        try:
            return pickle.loads(self.data[name])
        except:
            raise AttributeError

    @property
    def sample_uid(self):
        return self.experiment_sample.uid

    def as_dict(self):
        """Returns a flattened dictionary of the unpickled element values in
        VarantEvidence.data.
        """
        return get_flattened_unpickled_data(self.data)

    @classmethod
    def default_view_fields(clazz):
        return [
            {'field':'gt_type'},
            {'field':'sample_uid'},
        ]


###############################################################################
# Analysis
###############################################################################

class VariantToVariantSet(Model):
    """Relationship between variants and variant sets.

    In addition to linking a variant to a set, this model also allows
    strengthening the information content of the relationship by indicating
    which specific ExperimentSamples this relationship is valid for.
    """

    variant = models.ForeignKey('Variant')

    variant_set = models.ForeignKey('VariantSet')

    sample_variant_set_association = models.ManyToManyField('ExperimentSample',
            blank=True, null=True)

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'variant'},
                {'field':'sample_variant_set_association'}]


class VariantSet(Model):
    """Model for grouping together variants for analysis.

    This object can also be thought of a 'tag' for a set of variants. For
    example, we might create a VariantSet called 'c321D Designed Changes'
    to represent the set of Variants that were intended for mutation.

    Variants hold a list of Variant objects, and each can, but do not have to,
    point to one or more VariantToExperimentSample objects.

    Each variant set can contain variants from multiple alignments or samples,
    but all variants must belong to a single reference genome.

    TODO: In the future, we might come up with a framework for transferring
    variants or variant sets to new reference genomes via LiftOver or something
    similar.

    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(VariantSet)))

    label = models.CharField(max_length=256)

    reference_genome = models.ForeignKey('ReferenceGenome')

    variants = models.ManyToManyField('Variant', blank=True, null=True,
        # TODO: find correct syntax for limit_choices_to here
        #limit_choices_to = {'reference_genome' : self.reference_genome},
        through = 'VariantToVariantSet')

    # Datasets pointing to files on the system
    # Primarily for VCF files uploaded by the user to describe putative vars
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
        verbose_name="Datasets")

    def __unicode__(self):
        return self.label

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'label'},
                {'field':'reference_genome'}]

    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'genome_designer.main.views.variant_set_view',
                args=(self.reference_genome.project.uid, self.uid))

    def get_model_data_root(self):
        """Get the root location for all data of this type in the project.
        """
        return os.path.join(
                self.reference_genome.project.get_model_data_dir(),
                'variant_sets')

    def get_model_data_dir(self):
        """Get the full path to the location of this model's data.
        """
        return os.path.join(self.get_model_data_root(), str(self.uid))

    def ensure_model_data_dir_exists(self):
        """Ensure that a data directory exists for this model.
        """
        # Make sure the root exists.
        ensure_exists_0775_dir(self.get_model_data_root())

        # Check whether the data dir exists, and create it if not.
        return ensure_exists_0775_dir(self.get_model_data_dir())

# When a new Project is created, create the data directory.
def post_variant_set_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_variant_set_create, sender=VariantSet,
        dispatch_uid='variant_set_create')

class VariantFilter(Model):
    """Model used to save a string representation of a filter used to sift
    through Variants.
    """
    pass


class Region(Model):
    """Semantic annotation for a disjoint set of intervals in a
    ReferenceGenome.

    This allows the user to ask semantically deeper questions.
    """
    uid = models.CharField(max_length=36, unique=True,
            default=(lambda: short_uuid(Region)))

    reference_genome = models.ForeignKey('ReferenceGenome')

    # Human-readable identifier.
    # Depending on the type and how disciplined we are with development,
    # this could further be semantically meaningful (e.g. gene name).
    label = models.CharField('Sample Name', max_length=256)

    class TYPE:
        """The type of this dataset.

        Limit to 40-chars as per Dataset.type field def.
        """
        CALLABLE = 'c'
        GENE = 'g'
    TYPE_CHOICES = make_choices_tuple(TYPE)
    type = models.CharField(max_length=40, choices=TYPE_CHOICES)


class RegionInterval(Model):
    """One of possibly several intervals that describe a single region.
    """
    region = models.ForeignKey('Region')

    # One-indexed.
    start = models.BigIntegerField()

    # One-indexed.
    end = models.BigIntegerField()


class S3File(Model):
    """Model for keeping track of all files in S3 bucket.
    """
    bucket = models.CharField(max_length=200)

    # key is the actually name of the file stored in S3 bucket.
    key = models.CharField(max_length=200)

    # name is the original name of the file on uploader's machine
    name = models.CharField(max_length=200, null=True)
    created_at = models.DateTimeField(auto_now_add = True)

    def url(self):
        return "s3://%s/%s" % (self.bucket, self.key)

    def __unicode__(self):
        return unicode(self.url())
