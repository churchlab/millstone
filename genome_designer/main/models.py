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
"""
import os
import re
import stat
from uuid import uuid4

from django.contrib.auth.models import User
from django.db import models
from django.db.models import Model
from django.db.models.signals import post_save

import settings


###############################################################################
# Utilities
###############################################################################

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
        os.mkdir(destination)
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
    uid = models.CharField(max_length=36,
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
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(Dataset)))

    # The type of data this represents (e.g. Dataset.Type.BWA_ALIGN).
    # This is a semantic identifier for the kinds of operations
    # that can be performed with this Dataset.
    class TYPE:
        """The type of this dataset.

        Limit to 40-chars as per Dataset.type field def.
        """
        REFERENCE_GENOME_FASTA = 'rgf' # fasta
        REFERENCE_GENOME_GENBANK = 'rgg'
        FASTQ1 = 'f1'
        FASTQ2 = 'f2'
    TYPE_CHOICES = make_choices_tuple(TYPE)
    type = models.CharField(max_length=40, choices=TYPE_CHOICES)

    # Human-readable identifier. Also used for JBrowse.
    label = models.CharField(max_length=256)

    # Location on the filesystem relative to settings.MEDIA_ROOT.
    filesystem_location = models.CharField(max_length=512)

    def __unicode__(self):
        return self.label

    def get_absolute_location(self):
        """Returns the full path to the file on the filesystem.
        """
        return os.path.join(settings.MEDIA_ROOT, self.filesystem_location)

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
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(Project)))

    # The project owner.
    # TODO: Implement permissions system so that projects can be shared.
    owner = models.ForeignKey('UserProfile')

    # The human-readable title of the project.
    title = models.CharField(max_length=256)

    def __unicode__(self):
        return self.title + '-' + str(self.owner)

    def get_model_data_root(self):
        """Get the root location where all user data is stores."""
        return os.path.join(settings.MEDIA_ROOT, 'projects')

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

# When a new Project is created, create the data directory.
def post_project_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_project_create, sender=Project,
        dispatch_uid='project_create')


class ReferenceGenome(Model):
    """A reference genome relative to which alignments and analysis are
    performed.
    """
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(ReferenceGenome)))

    # A ReferenceGenome belongs to a single Project.
    project = models.ForeignKey('Project')

    # A human-readable label for this genome.
    label = models.CharField(max_length=256)

    # Number of chromosomes.
    num_chromosomes = models.IntegerField()

    # Number of chromosomes.
    num_bases = models.BigIntegerField()

    # Datasets pointing to files on the system (e.g. .fasta files, etc.)
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True)

    def __unicode__(self):
        return self.label

    def get_model_data_root(self):
        """Get the root location where all user data is stores."""
        return os.path.join(self.project.get_model_data_dir(), 'ref_genomes')

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

# When a new ReferenceGenome is created, create its data dir.
def post_ref_genome_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_ref_genome_create, sender=ReferenceGenome,
        dispatch_uid='ref_genome_create')


class ExperimentSample(Model):
    """Model representing data for a particular experiment sample.

    Usually this corresponds to a pair of fastq reads for a particular colony,
    after de-multiplexing.
    """
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(ExperimentSample)))

    # A Sample belongs to a single Project.
    project = models.ForeignKey('Project')

    # Human-readable identifier.
    label = models.CharField(max_length=256, blank=True)

    # The datasets associated with this sample. The semantic sense of the
    # dataset can be determined from the Dataset.type field.
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True)


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
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(AlignmentGroup)))

    # Human-readable identifier.
    label = models.CharField(max_length=256, blank=True)

    # All alignments in this set are relative to this genome.
    reference_genome = models.ForeignKey('ReferenceGenome')

    # The aligner tool used for this alignment.
    class ALIGNER:
        """Constants for representing the aligner type.
        """
        BOWTIE2 = 'BOWTIE2'
        BWA = 'BWA'
    ALIGNER_CHOICES = make_choices_tuple(ALIGNER)
    aligner = models.CharField(max_length=10, choices=ALIGNER_CHOICES)

    # Times for the alignment run.
    start_time = models.DateTimeField(auto_now=True)
    end_time = models.DateTimeField(auto_now=True)


class ExperimentSampleToAlignment(Model):
    """Model that describes the alignment of a single ExperimentSample
    to an Alignment.
    """
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(ExperimentSampleToAlignment)))

    alignment_group = models.ForeignKey('AlignmentGroup')

    experiment_sample = models.ForeignKey('ExperimentSample')


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
    uid = models.CharField(max_length=36,
            default=(lambda: short_uuid(Variant)))

    class TYPE:
        DELETION = 'DELETION'
        INSERTION = 'INSERTION'
        TRANSITION = 'TRANSITION'
        TRANSVERSION = 'TRANSVERSION'
        COMPLEX = 'COMPLEX' # Multi-base in different genomes
    TYPE_CHOICES = make_choices_tuple(TYPE)
    type = models.CharField(max_length=40, choices=TYPE_CHOICES)

    reference_genome = models.ForeignKey('ReferenceGenome')

    chromosome = models.CharField(max_length=256, blank=True)

    position = models.BigIntegerField()

    ref_value = models.TextField();

    alt_value = models.TextField();


class VariantToExperimentSample(Model):
    """Model that represents the relationship between a particular Variant
    and a particular ExperimentSample.
    """
    pass


class VariantCallerCommonData(Model):
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
    pass


class VariantEvidence(Model):
    """Evidence for a particular variant occurring in a particular
    ExperimentSample.
    """
    pass


###############################################################################
# Analysis
###############################################################################

class VariantSet(Model):
    """Model for grouping together variants for analysis.

    This object can also be thought of a 'tag' for a set of variants. For
    example, we might create a VariantSet called 'c321D Designed Changes'
    to represent the set of Variants that were intended for mutation.
    """
    pass


class VariantFilter(Model):
    """Model used to save a string representation of a filter used to sift
    through Variants.
    """
    pass
