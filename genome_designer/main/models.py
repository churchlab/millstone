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

from contextlib import contextmanager
import os
import shutil
import subprocess

from custom_fields import PostgresJsonField
from django.conf import settings
from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.db import models
from django.db.models import Model

from model_utils import assert_unique_types
from model_utils import ensure_exists_0775_dir
from model_utils import get_dataset_with_type
from model_utils import make_choices_tuple
from model_utils import UniqueUidModelMixin
from model_utils import VisibleFieldMixin
from settings import TOOLS_DIR
from scripts.filter_key_map_constants import MAP_KEY__ALTERNATE
from scripts.filter_key_map_constants import MAP_KEY__COMMON_DATA
from scripts.filter_key_map_constants import MAP_KEY__EVIDENCE
from scripts.filter_key_map_constants import MAP_KEY__EXPERIMENT_SAMPLE
from scripts.util import uppercase_underscore

BGZIP_BINARY = '%s/tabix/bgzip' % TOOLS_DIR


###############################################################################
# User-related models
###############################################################################

class UserProfile(UniqueUidModelMixin):
    """A UserProfile which is separate from the django auth User.

    This references the auth.User and opens up the possibility of
    adding additional fields.
    """
    # A one-to-one mapping to the django User model.
    user = models.OneToOneField(User)

    def __unicode__(self):
        return self.user.username


###############################################################################
# Data wrappers
###############################################################################

class Dataset(UniqueUidModelMixin):
    """A specific data file with a location on the filesystem.

    Basically a wrapper for a file on the file system.

    This is similar to the Galaxy notion of a dataset.
    """

    # The type of data this represents (e.g. Dataset.Type.BWA_ALIGN).
    # This is a semantic identifier for the kinds of operations
    # that can be performed with this Dataset.
    class TYPE:
        """The type of this dataset.

        Limit to 40-chars as per Dataset.type field def.

        For internal strings, we will convert to ALL_CAPS_W_UNDERSCORES.

        """
        REFERENCE_GENOME_FASTA = 'Reference Genome FASTA' # fasta
        REFERENCE_GENOME_GENBANK = 'Reference Genome Genbank' #genbank
        REFERENCE_GENOME_GFF = 'Reference Genome GFF' #gff
        FASTQ1 = 'FASTQ Forward'
        FASTQ2 = 'FASTQ Reverse'
        BWA_ALIGN = 'BWA BAM'
        BWA_ALIGN_ERROR = 'BWA Alignment Error'
        VCF_FREEBAYES = 'Freebayes VCF'
        VCF_PINDEL = 'Pindel VCF'
        VCF_DELLY = 'Delly VCF'
        VCF_USERINPUT = 'User VCF'
        VCF_FREEBAYES_SNPEFF = 'SNPEff VCF'
        BED_CALLABLE_LOCI = 'Flagged Regions BED'
        PICARD_INSERT_METRICS = 'Picard Insert'


    TYPE_CHOICES = make_choices_tuple(TYPE)
    type = models.CharField(max_length=40, choices=TYPE_CHOICES)

    # This relationship lets us know where the dataset points. This
    # is important in case we want to duplicate this dataset in order
    # to make a compressed/uncompressed version - we need to hook it
    # up to the correct related models after copying. 
    TYPE_TO_RELATED_MODEL = {
        TYPE.REFERENCE_GENOME_FASTA : 'referencegenome_set',
        TYPE.REFERENCE_GENOME_GENBANK : 'referencegenome_set',
        TYPE.FASTQ1 : 'experimentsample_set',
        TYPE.FASTQ2 : 'experimentsample_set',
        TYPE.BWA_ALIGN : 'experimentsampletoalignment_set',
        TYPE.BWA_ALIGN_ERROR : 'alignmentgroup_set',
        TYPE.VCF_FREEBAYES : 'alignmentgroup_set',
        TYPE.VCF_PINDEL : 'alignmentgroup_set',
        TYPE.VCF_DELLY : 'alignmentgroup_set',
        TYPE.VCF_USERINPUT : 'variantset_set',
        TYPE.VCF_FREEBAYES_SNPEFF : 'alignmentgroup_set'
    }

    # Human-readable identifier. Also used for JBrowse.
    label = models.CharField(max_length=256)

    # Location on the filesystem relative to settings.MEDIA_ROOT.
    filesystem_location = models.CharField(max_length=512, blank=True)

    # Associated with a separate index? (e.g. for vcf/tabix and bam files)
    filesystem_idx_location = models.CharField(max_length=512, blank=True)

    # When the dataset is a result of a computation, we'll set a status on it.
    # NOTE: The reliability of the present implementation of this model feature
    # is questionable.
    class STATUS:
        """
        The status of running this Dataset.

        Limit to 40-chars as per Dataset.status field def.
        """
        UNKNOWN = 'UNKNOWN'
        NOT_STARTED = 'NOT_STARTED'
        COMPUTING = 'COMPUTING'
        ALIGNING = 'ALIGNING'
        VARIANT_CALLING = 'VARIANT_CALLING'
        READY = 'READY'
        FAILED = 'FAILED'
        COPYING = 'COPYING'
        QUEUED_TO_COPY = 'QUEUED_TO_COPY'
        VERIFYING = 'VERIFYING'
    STATUS_CHOICES = make_choices_tuple(STATUS)
    status = models.CharField(max_length=40, choices=STATUS_CHOICES,
            default=STATUS.READY)

    # Dictionary of compression suffixes and programs to use to perform
    # various actions on a pipe
    COMPRESSION_TYPES = {
        '.gz': {'cat': ('gzip','-dc'), 'zip': ('gzip','-c')},
        '.bz2': {'cat': ('bzcat',), 'zip': ('bzip2','-c')},
        '.zip': {'cat': ('unzip','-p'), 'zip': ('zip','-')},
        '.bgz': {'cat': (BGZIP_BINARY,'-dc'), 'zip': (BGZIP_BINARY,'-c')},
    }

    def __unicode__(self):
        return self.label

    def get_absolute_location(self):
        """Returns the full path to the file on the filesystem.
        """
        return os.path.join(settings.PWD, settings.MEDIA_ROOT,
                self.filesystem_location)

    def get_absolute_idx_location(self):
        return os.path.join(settings.PWD, settings.MEDIA_ROOT,
                self.filesystem_idx_location)

    def is_compressed(self):
        """
        Checks file path for .bz2 or .gz ending, and if so, returns true.
        """
        return self.filesystem_location.endswith(
                tuple(self.COMPRESSION_TYPES.keys()))

    def is_indexed(self):
        """
        Checks if dataset has idx location.
        """
        return not (self.filesystem_idx_location == '')

    def wrap_if_compressed(self):
        """ This helper function returns a process substitution string
        to be used by /bin/bash if the fastq read file is compressed, otherwise
        it just returns get_absolute_location().
        """
        absolute_location = self.get_absolute_location()
        if self.is_compressed():
            extension = os.path.splitext(self.filesystem_location)[1]
            program = ' '.join(self.COMPRESSION_TYPES[extension]['cat'])
            return '<({:s} {:s})'.format(
                    program, absolute_location)
        else:
            return absolute_location

    def internal_string(self, parent_entity):
        """
        A string used internally to describe a dataset for an entity.
        Our convention is 
            parent_entity.uid + '_' + dataset.type 
            (uppercased, whitespace as underscores)
        """
        return str(parent_entity.uid) + '_' + uppercase_underscore(self.type)

    def external_string(self, parent_entity):
        """
        A string used externally to describe a dataset for an entity.
        Our convention is 
            parent_entity.label + ' ' + dataset.type
        """
        return str(parent_entity.label) + ' ' + self.type


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

    def get_related_model_set(self):
        return getattr(self, Dataset.TYPE_TO_RELATED_MODEL[self.type])

    def make_compressed(self, compression_type):
        """
        Generate a new compressed version of this dataset.

        For some cases (like generating a compressed TABIX-indexed VCF),
        we want to take a dataset and generate a compressed version of
        the file (as a separate model instance) with the same associations 
        to other related model instances. 

        TODO: We could just replace the uncompressed version with the
        compressed version with the compressed version, but right now that's
        too much work, since we'd need to check every time to see if the file
        was compressed, and depending on the tool, decide if we'd need to
        decompress it via pipe, or write the decompressed version as a new
        file, etc etc. So, currently the replace option is unimplemented.
        """
        # Check that compression_type is valid
        assert compression_type in Dataset.COMPRESSION_TYPES, (
                'compression_type is invalid, {:s} is not one of: {:s}'.format(
                        compression_type, Dataset.COMPRESSION_TYPES.keys()))

        # Check that this dataset isn't itself already compressed
        assert self.is_compressed() is False, (
                'This dataset is already compressed.')

        # Check that a compressed dataset isn't already associated with a 
        # related model (probably just one related model).    
        related_models = self.get_related_model_set().all()
        for obj in related_models:
            old_compressed_dataset = get_dataset_with_type(obj, self.type, 
                    compressed=True)
            # TODO: In this case, do we want to just return this?
            # Maybe with a warning?
            assert old_compressed_dataset is None, (
                'A related model already compressed' +
                'this dataset: {:s}'.format(
                        compressed_dataset.filesystem_location))
            

        # Generate the new compressed dataset file 
        # by adding the compression_type suffix
        orig_file = self.get_absolute_location()
        new_compressed_file = orig_file + compression_type

        compression_command = Dataset.COMPRESSION_TYPES[
                compression_type]['zip']

        with open(orig_file, 'rb') as fh_in:
            with open(new_compressed_file, 'wb') as fh_out:
                subprocess.check_call(compression_command, 
                        stdin=fh_in, stdout=fh_out)

        # Generate the new dataset model object
        # need relative path, not absolute
        new_compressed_file_rel = self.filesystem_location + compression_type

        new_compressed_dataset = Dataset.objects.create(
            label= self.label + ' (compressed)', 
            type= self.type, 
            filesystem_location= new_compressed_file_rel)

        # Finally, add this new compressed dataset to the dataset_set
        # field in all the related model objects that point to the
        # uncompressed version
        [obj.dataset_set.add(new_compressed_dataset) for obj in related_models]

        return new_compressed_dataset


# Make sure the Dataset types are unique. This runs once at startup.
assert_unique_types(Dataset.TYPE)

###############################################################################
# Project models
###############################################################################

class Project(UniqueUidModelMixin):
    """A single project belonging to a user.

    A project groups together ReferenceGenomes, ExperimentSamples, and other
    data generated by tools during analysis.
    """
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


class ReferenceGenome(UniqueUidModelMixin):
    """A reference genome relative to which alignments and analysis are
    performed.
    """

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

    # a key/value list of all possible VCF and sample metadata fields, stored
    # as a JsonField and dynamically updated by dynamic_snp_filter_key_map.py
    variant_key_map = PostgresJsonField()

    # Bit that indicates whether the materialized view is up to date.
    # This design decision puts a lot on the developer to remember to set this
    # false whenever any data changes that would require a refresh of the
    # materialized view.
    is_materialized_variant_view_valid = models.BooleanField(default=False)

    def __unicode__(self):
        return self.label

    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'main.views.reference_genome_view',
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
                    maybe_force_nginx + '/jbrowse/gd_data/',
                    'projects',
                    str(self.project.uid),
                    'ref_genomes',
                    str(self.uid),
                    'jbrowse')

    def get_client_jbrowse_link(self):
        """Returns the link to jbrowse redirect for this ReferenceGenome.

        Example url for user with uid 'abc', and project id 'xyz', and
        refgenome id 456:
            '/redirect_jbrowse?data=gd_data/abc/projects/xyz/ref_genomes/456/jbrowse/'
        """
        return '/redirect_jbrowse?data=' + self.get_client_jbrowse_data_path()

    def is_annotated(self):
        """For several steps (notably snpEff), we want to check that this
        ReferenceGenome is annotated (i.e. it has a genbank file associated
        with it). This function returns True if a genbank file is available.
        """
        return self.dataset_set.filter(
                type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).exists()

    def get_variant_caller_common_map(self):
        return self.variant_key_map[MAP_KEY__COMMON_DATA]

    def get_variant_alternate_map(self):
        return self.variant_key_map[MAP_KEY__ALTERNATE]

    def get_variant_evidence_map(self):
        return self.variant_key_map[MAP_KEY__EVIDENCE]

    def get_experiment_sample_map(self):
        return self.variant_key_map[MAP_KEY__EXPERIMENT_SAMPLE]


    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [{'field':'label'},
                {'field':'num_chromosomes', 'verbose':'# Chromosomes'},
                {'field':'num_bases', 'verbose':'Total Size'}]

    def invalidate_materialized_view(self):
        self.is_materialized_variant_view_valid = False
        self.save(update_fields=['is_materialized_variant_view_valid'])


class ExperimentSample(UniqueUidModelMixin):
    """Model representing data for a particular experiment sample.

    Usually this corresponds to a pair of fastq reads for a particular
    bacterial clone or colony, after barcode removal/de-multiplexing.
    """

    # A Sample belongs to a single Project.
    project = models.ForeignKey('Project')

    # Human-readable identifier.
    label = models.CharField('Sample Name', max_length=256)

    # The datasets associated with this sample. The semantic sense of the
    # dataset can be determined from the Dataset.type field.
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
        verbose_name="Datasets")

    # User speciified data fields corresponding to the sample.
    # Examples: Growth rate, GFP amount, phenotype, # of mage cycles, etc.
    data = PostgresJsonField()

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
            return self.data[name]
        except:
            raise AttributeError

    @property
    def status(self):
        """The status of the data underlying this data.
        """
        status_string = 'NO_DATA'
        fastq1_dataset_queryset = self.dataset_set.filter(
                type=Dataset.TYPE.FASTQ1)
        if len(fastq1_dataset_queryset) > 1:
            return 'ERROR: More than one forward reads source'
        if len(fastq1_dataset_queryset) == 1:
            status_string = 'FASTQ1: %s' % fastq1_dataset_queryset[0].status
            # Maybe add reverse reads.
            fastq2_dataset_queryset = self.dataset_set.filter(
                    type=Dataset.TYPE.FASTQ2)
            if len(fastq2_dataset_queryset) > 1:
                return 'ERROR: More than one reverse reads source'
            if len(fastq2_dataset_queryset) == 1:
                status_string += (
                        ' | FASTQ2:  %s' % fastq2_dataset_queryset[0].status)
        return status_string

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

        return [
            {'field': 'label'},
            {'field': 'status'},
            {'field': 'uid', 'verbose':'Internal ID'}]


class AlignmentGroup(UniqueUidModelMixin):
    """Collection of alignments of several related ExperimentSamples to the
    same ReferenceGenome.

    The reason for grouping alignments together is that our variant operations
    are generally relative to a single reference genome, and further, variant
    calling tools often take multiple alignments as input, thus it makes sense
    to group these in the database.

    For a one-to-one mapping of Alignment to Sample, see
    ExperimentSampleToAlignment.
    """

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
    start_time = models.DateTimeField(blank=True, null=True)
    end_time = models.DateTimeField(blank=True, null=True)

    # Datasets pointing to files on the system (e.g. .fasta files, etc.)
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
            verbose_name="Datasets")

    class STATUS:
        """
        The status of running this Dataset.

        Limit to 40-chars as per Dataset.status field def.
        """
        NOT_STARTED = 'NOT_STARTED'
        ALIGNING = 'ALIGNING'
        VARIANT_CALLING = 'VARIANT_CALLING'
        COMPLETED = 'COMPLETED'
        FAILED = 'FAILED'
        UNKNOWN = 'UNKNOWN'
    STATUS_CHOICES = make_choices_tuple(STATUS)
    
    status = models.CharField('Alignment Status', 
            max_length=40, choices=STATUS_CHOICES, default=STATUS.NOT_STARTED)

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
                'main.views.alignment_view',
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

    def get_samples(self):
        """Many different tasks require getting the sample (or their UIDs)
        that are in this alignment group.
        """
        return ExperimentSample.objects.filter(
                experimentsampletoalignment__alignment_group=self)


class ExperimentSampleToAlignment(UniqueUidModelMixin):
    """Model that describes the alignment of a single ExperimentSample
    to an AlignmentGroup.
    """

    alignment_group = models.ForeignKey('AlignmentGroup')

    experiment_sample = models.ForeignKey('ExperimentSample')

    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True)

    @property
    def status(self):
        """The status of a running alignment job.
        """
        alignment_datasets = self.dataset_set.filter(
                type=Dataset.TYPE.BWA_ALIGN)
        assert len(alignment_datasets) <= 1, (
                "Expected only one alignment dataset.")
        if len(alignment_datasets) == 1:
            return alignment_datasets[0].status
        return 'UNDEFINED'

    @property
    def error_link(self):
        return ('<a href="' +
                reverse(
                        'main.views.sample_alignment_error_view',
                        args=(self.alignment_group.reference_genome.project.uid,
                                self.alignment_group.uid,
                                self.uid)) +
                        '">log output</a>')

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [
            {'field': 'experiment_sample'},
            {'field': 'status', 'verbose': 'Job Status'},
            {'field': 'error_link', 'verbose': 'Sample Alignment Log', 'is_href': True},
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


###############################################################################
# Variants (SNVs and SVs)
###############################################################################

class Variant(UniqueUidModelMixin):
    """An instance of a variation relative to a reference genome.

    This might be, for example, a SNV (single nucleotide variation) or a bigger
    SV (structural variation). We are intentionally using a unified model for
    these two classes of variations as the fundamental unit of genome analysis
    is really a diff.

    TODO: See code from Gemini paper (PLOS ONE 7/18/13) for ideas.

    A variant need not necessarily be associated with a specific sample; the
    VariantToExperimentSample model handles this association.
    """

    class TYPE:
        DELETION = 'DELETION'
        INSERTION = 'INSERTION'
        TRANSITION = 'TRANSITION'
        TRANSVERSION = 'TRANSVERSION'
        DUPLICATION = 'DUPLICATION'
        INVERSION = 'INVERSION'
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
                'main.views.single_variant_view',
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
    data = PostgresJsonField()

    alignment_group = models.ForeignKey('AlignmentGroup')

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
            return self.data[name]
        except:
            raise AttributeError

    @classmethod
    def default_view_fields(clazz):
        return []


class VariantAlternate(UniqueUidModelMixin, VisibleFieldMixin):
    """A model listing alternate alleles for each variant."""

    # Null is true here because we are adding this relationship during Variant's
    # overloaded __init__() so it hasn't been saved() yet. Otherwise it throws
    # an django.db.utils.IntegrityError: 
    #    main_variantalternate.variant_id may not be NULL
    variant = models.ForeignKey('Variant', null=True)

    alt_value = models.TextField('Alt')

    is_primary = models.BooleanField(default='False')

    # this json fields holds all PER ALT data (INFO data with num -1)
    data = PostgresJsonField()

    def __unicode__(self):
        return 'var: ' + str(self.variant) + ', alt:' + self.alt_value

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


class VariantEvidence(UniqueUidModelMixin, VisibleFieldMixin):
    """
    Evidence for a particular variant occurring in a particular
    ExperimentSample.
    """

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
    data = PostgresJsonField()

    def __init__(self, *args, **kwargs):
        # HACK: Manually cache data to avoid expensive lookups.
        self.manually_cached_data = {}

        # call super's __init__ without the alt_value field if present
        super(VariantEvidence, self).__init__(*args, **kwargs)


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
            return self.data[name]
        except:
            raise AttributeError


    def create_variant_alternate_association(self):
        gt_bases = self.data['GT_BASES']

        # If this variant evidence is a non-call, no need to add alt alleles.
        if gt_bases is None:
            return

        assert ('|' not in gt_bases), (
                'GT bases string is phased;' +
                'this is not handled and should never happen...')

        # The gt_bases string looks like, e.g. 'A/AT'. Loop over alts.
        for gt_base in gt_bases.split('/'):
            try:
                variant = self.variant_caller_common_data.variant

                # Skip if this is not an alternate allele
                if variant.ref_value == gt_base:
                    continue

                self.variantalternate_set.add(
                        VariantAlternate.objects.get(
                            variant=variant,
                            alt_value=gt_base
                ))

            except VariantAlternate.DoesNotExist:
                # Should not happen.
                print ('Attempt to add a SampleEvidence with an alternate ' +
                        'allele that is not present for this variant!')
                raise

    @property
    def sample_uid(self):
        if 'sample_uid' in self.manually_cached_data:
            return self.manually_cached_data['sample_uid']

        # Otherwise, probably do DB lookup. Guarantee correctness.
        return self.experiment_sample.uid

    @classmethod
    def default_view_fields(clazz):
        return [
            {'field':'gt_type'},
            {'field':'sample_uid', 'verbose':'Samples'},
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


class VariantSet(UniqueUidModelMixin):
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
                'main.views.variant_set_view',
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


class Region(UniqueUidModelMixin):
    """Semantic annotation for a disjoint set of intervals in a
    ReferenceGenome.

    This allows the user to ask semantically deeper questions.
    """
    reference_genome = models.ForeignKey('ReferenceGenome')

    # Human-readable identifier.
    # Depending on the type and how disciplined we are with development,
    # this could further be semantically meaningful (e.g. gene name).
    label = models.CharField('Region Name', max_length=256)

    class TYPE:
        """The type of this region.

        Limit to 40-chars as per Dataset.type field def.
        """
        POOR_MAPPING_QUALITY = 'POOR_MAPPING_QUALITY'
        GENE = 'GENE'

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
