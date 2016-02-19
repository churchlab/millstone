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
from datetime import datetime
import json
import os
import re
import shutil
import subprocess

from custom_fields import PostgresJsonField
from django.conf import settings
from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.db import models
from django.db.models import Model

# from genome_finish.insertion_placement_read_trkg import Junction
from genome_finish.contig_display_utils import create_contig_junction_links
from model_utils import assert_unique_types
from model_utils import ensure_exists_0775_dir
from model_utils import get_dataset_with_type
from model_utils import make_choices_tuple
from model_utils import UniqueUidModelMixin
from model_utils import VisibleFieldMixin
from utils import uppercase_underscore
from variants.filter_key_map_constants import MAP_KEY__ALTERNATE
from variants.filter_key_map_constants import MAP_KEY__COMMON_DATA
from variants.filter_key_map_constants import MAP_KEY__EVIDENCE
from variants.filter_key_map_constants import MAP_KEY__EXPERIMENT_SAMPLE
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager


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
        BWA_ALTALIGN = 'BWA Alternate Alignment Reads'
        BWA_DISCORDANT = 'BWA BAM Discordant Paired Reads'
        BWA_SPLIT = 'BWA BAM Split Reads'
        BWA_UNMAPPED = 'BWA Unmapped Reads'
        BWA_CLIPPED = 'BWA Clipped Reads'
        BWA_PILED = 'BWA Piled Reads'
        BWA_SV_INDICANTS = 'BWA Structural Variant Indicating Reads'
        BWA_FOR_DE_NOVO_ASSEMBLY = 'BWA for De Novo Assembly'
        BWA_ALIGN_ERROR = 'BWA Alignment Error'
        VCF_FREEBAYES = 'Freebayes VCF'
        VCF_PINDEL = 'Pindel VCF'
        VCF_DELLY = 'Delly VCF'
        VCF_LUMPY = 'Lumpy VCF'
        VCF_USERINPUT = 'User VCF'
        VCF_FREEBAYES_SNPEFF = 'SNPEff VCF'
        VCF_LUMPY_SNPEFF = 'Lumpy SNPEff VCF'
        VCF_PINDEL_SNPEFF = 'Pindel SNPEff VCF'
        VCF_DE_NOVO_ASSEMBLED_CONTIGS = 'De Novo Assembled Contigs VCF'
        VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK = 'De Novo Assembly Graph Walk VCF'
        VCF_COV_DETECT_DELETIONS = 'Deletions Detected by Coverage Evidence'
        BED_CALLABLE_LOCI = 'Flagged Regions BED'
        LUMPY_INSERT_METRICS_HISTOGRAM = 'Lumpy Insert Metrics Histogram'
        LUMPY_INSERT_METRICS_MEAN_STDEV = 'Lumpy Insert Metrics Mean Stdev'
        FASTQC1_HTML = 'FASTQC Forward HTML Output'
        FASTQC2_HTML = 'FASTQC Reverse HTML Output'

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
        TYPE.BWA_ALTALIGN : 'experimentsampletoalignment_set',
        TYPE.BWA_DISCORDANT : 'experimentsampletoalignment_set',
        TYPE.BWA_SPLIT : 'experimentsampletoalignment_set',
        TYPE.BWA_CLIPPED : 'experimentsampletoalignment_set',
        TYPE.BWA_PILED : 'experimentsampletoalignment_set',
        TYPE.BWA_SV_INDICANTS : 'experimentsampletoalignment_set',
        TYPE.BWA_ALIGN_ERROR : 'alignmentgroup_set',
        TYPE.VCF_FREEBAYES : 'alignmentgroup_set',
        TYPE.VCF_PINDEL : 'alignmentgroup_set',
        TYPE.VCF_LUMPY : 'alignmentgroup_set',
        TYPE.VCF_DELLY : 'alignmentgroup_set',
        TYPE.VCF_LUMPY : 'alignmentgroup_set',
        TYPE.VCF_LUMPY_SNPEFF: 'alignmentgroup_set',
        TYPE.VCF_USERINPUT : 'variantset_set',
        TYPE.VCF_FREEBAYES_SNPEFF : 'alignmentgroup_set',
        TYPE.VCF_DE_NOVO_ASSEMBLED_CONTIGS : 'experimentsampletoalignment_set',
        TYPE.VCF_DE_NOVO_ASSEMBLY_GRAPH_WALK : (
                'experimentsampletoalignment_set'),
        TYPE.FASTQC1_HTML: 'experimentsample_set',
        TYPE.FASTQC2_HTML: 'experimentsample_set',
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
        QC = 'RUNNING_QC'
        AWAITING_UPLOAD = 'AWAITING_UPLOAD'
    STATUS_CHOICES = make_choices_tuple(STATUS)
    status = models.CharField(max_length=40, choices=STATUS_CHOICES,
            default=STATUS.READY)

    # Dictionary of compression suffixes and programs to use to perform
    # various actions on a pipe
    COMPRESSION_TYPES = {
        '.gz': {
            'cat': ('gzip', '-dc'),
            'zip': ('gzip', '-c')
        },
        '.bz2': {
            'cat': ('bzcat',),
            'zip': ('bzip2', '-c')
        },
        '.zip': {
            'cat': ('unzip', '-p'),
            'zip': ('zip', '-')
        },
        '.bgz': {
            'cat': (settings.BGZIP_BINARY, '-dc'),
            'zip': (settings.BGZIP_BINARY, '-c')
        },
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

    def delete_underlying_data(self):
        """Deletes data from filesystem.
        """
        full_fs_location = self.get_absolute_location()
        if os.path.exists(full_fs_location):
            os.remove(full_fs_location)

        full_fs_index_location = self.get_absolute_idx_location()
        if os.path.exists(full_fs_index_location):
            os.remove(full_fs_index_location)

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


class Chromosome(UniqueUidModelMixin):
    """A locus belonging to a reference genome which Variants 
    hold foreign keys to.  May be a literal chromosome,
    bacterial genome, or plasmid.
    """

    # A chromosome belongs to a single ReferenceGenome
    reference_genome = models.ForeignKey('ReferenceGenome')

    # Chromosome label
    label = models.CharField(verbose_name="Name", max_length=256)

    # The unique id of the SeqRecord object corresponding to this Chromosome.
    # In a genbank/multi-fasta file, the sequence belonging to each chromosome
    # carries a unique identifier which is parsed by BioPython's SeqIO module
    # as the .id attribute of the SeqRecord object.  This field ties our
    # Chromosome model to a specific chromosome in a reference genome
    # fasta/genbank dataset.  The seqrecord_id field does not necesarilly
    # carry any comprehensible information about the Chromosome, it is only an
    # identifier, and the description of the Chromosome is given by the
    # label field.
    # Ex: A reporter plasmid Chromosome:
    #       seqrecord_id: pl1839
    #       label: Reporter plasmid carrying RFP on a lac promoter
    seqrecord_id = models.CharField(
            verbose_name="SeqRecord ID", max_length=256, default="chrom_1")

    # Number of bases on the Chromosome
    num_bases = models.BigIntegerField()

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [
            {'field': 'label', 'verbose': 'Chromosome Name'},
            {'field': 'num_bases', 'verbose':'Bases'},
            {'field': 'uid'}
        ]


class ReferenceGenome(UniqueUidModelMixin):
    """A reference genome relative to which alignments and analysis are
    performed.
    """

    # A ReferenceGenome belongs to a single Project.
    project = models.ForeignKey('Project')

    # A human-readable label for this genome.
    label = models.CharField(verbose_name="Name", max_length=256)

    # Datasets pointing to files on the system (e.g. .fasta files, etc.)
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
        verbose_name="Datasets")

    # a key/value list of all possible VCF and sample metadata fields, stored
    # as a JsonField and dynamically updated by dynamic_snp_filter_key_map.py
    variant_key_map = PostgresJsonField()

    # reference genome metadata field for storing key-value pairs of reference
    # genome related information e.g. metadata['is_from_de_novo_assembly']=True
    metadata = PostgresJsonField()

    # Bit that indicates whether the materialized view is up to date.
    # This design decision puts a lot on the developer to remember to set this
    # false whenever any data changes that would require a refresh of the
    # materialized view.
    is_materialized_variant_view_valid = models.BooleanField(default=False)

    def __unicode__(self):
        return self.label

    @property
    def num_chromosomes(self):
        """Number of Chromosomes belonging to the ReferenceGenome
        """
        return len(Chromosome.objects.filter(reference_genome = self))

    @property
    def num_bases(self):
        """Total number of bases of all Chromosomes belonging to
        the ReferenceGenome
        """
        return sum([chrom.num_bases for chrom in Chromosome.objects.filter(reference_genome = self)])

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
            assert False, "url is incorrect."
            return os.path.join(
                    'http://%s.s3.amazonaws.com/' % settings.S3_BUCKET,
                    'projects',
                    str(self.project.uid),
                    'ref_genomes',
                    str(self.uid),
                    'jbrowse')
        else:
            return os.path.join(
                    '/jbrowse/gd_data/',
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
            {'field': 'num_chromosomes', 'verbose': '# Chromosomes'},
            {'field': 'num_bases', 'verbose': 'Total Size'}
        ]

    def invalidate_materialized_view(self):
        self.is_materialized_variant_view_valid = False
        self.save(update_fields=['is_materialized_variant_view_valid'])

    def drop_materialized_view(self):
        """Deletes associated materialized view.
        """
        mvm = MeltedVariantMaterializedViewManager(self)
        mvm.drop()


class Contig(UniqueUidModelMixin):

    # A human-readable label for this genome.
    label = models.CharField(verbose_name="Name", max_length=256)

    # Number of bases in the Contig
    num_bases = models.BigIntegerField(default=0)

    # Datasets pointing to files on the system (e.g. .fasta files, etc.)
    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True,
        verbose_name="Datasets")

    # Reference genome which the insertion belongs to
    parent_reference_genome = models.ForeignKey('ReferenceGenome',
            related_name='+')

    # The sample alignment that provides evidence for the insertion
    experiment_sample_to_alignment = models.ForeignKey(
            'ExperimentSampleToAlignment')

    # The variant caller common data object associated
    variant_caller_common_data = models.ForeignKey('VariantCallerCommonData',
            blank=True, null=True)

    # Contig metadata field for storing key-value pairs of contig
    # related information e.g. metadata['is_from_de_novo_assembly']=True
    metadata = PostgresJsonField()

    def __getattr__(self, name):
        """Automatically called if an attribute is not found in the typical
        place.

        Our implementation checks the metadata dict, raises AttributeError if
        not found
        """
        try:
            return self.metadata[name]
        except:
            raise AttributeError

    def get_model_data_root(self):
        """Get the root location for all data of this type in the project.
        """
        return os.path.join(
                self.parent_reference_genome.project.get_model_data_dir(),
                'contigs')

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
        Contig.
        """
        return os.path.join(self.get_model_data_dir(), 'jbrowse')

    def ensure_jbrowse_dir(self):
        """Ensures that the jbrowse data dir exists."""
        return ensure_exists_0775_dir(self.get_jbrowse_directory_path())

    def get_client_jbrowse_data_path(self):
        if self.parent_reference_genome.project.is_s3_backed():
            assert False, "url is incorrect."
        else:
            return os.path.join(
                    '/jbrowse/gd_data/',
                    'projects',
                    str(self.parent_reference_genome.project.uid),
                    'contigs',
                    str(self.uid),
                    'jbrowse')

    def get_client_jbrowse_link(self):
        """Returns the link to jbrowse redirect for this Contig.
        """
        bam_dataset = self.dataset_set.get(type=Dataset.TYPE.BWA_ALIGN)
        bam_label = bam_dataset.internal_string(self)
        coverage_label = bam_dataset.internal_string(self) + '_COVERAGE'
        track_labels = (settings.JBROWSE_DEFAULT_TRACKS +
                [bam_label, coverage_label])

        link = '/redirect_jbrowse?data=' + self.get_client_jbrowse_data_path()
        link += '&tracks=' + ','.join(track_labels)
        return link

    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'main.views.contig_view',
                args=(self.parent_reference_genome.project.uid, self.uid))

    @property
    def coverage(self):
        return self.metadata.get('coverage', '')

    @property
    def chromosome(self):
        return self.metadata.get('chromosome', '')

    def get_contig_reads_track(self):
        bam_dataset = get_dataset_with_type(
                self,
                Dataset.TYPE.BWA_SV_INDICANTS)
        return str(bam_dataset.internal_string(self))

    @property
    def left_junctions_html(self):
        junctions = self.metadata.get('left_junctions', '')
        return create_contig_junction_links(self, junctions)

    @property
    def right_junctions_html(self):
        junctions = self.metadata.get('right_junctions', '')
        return create_contig_junction_links(self, junctions)

    @property
    def experiment_sample(self):
        return self.experiment_sample_to_alignment.experiment_sample.label

    @classmethod
    def get_field_order(clazz, **kwargs):
        """Get the order of the models for displaying on the front-end.
        Called by the adapter.
        """
        return [
            {'field': 'label'},
            {'field': 'experiment_sample'},
            {'field': 'num_bases', 'verbose': 'Contig Length'},
            {'field': 'coverage', 'verbose': 'Average Coverage'},
            {'field': 'chromosome'},
            {'field': 'left_junctions_html', 'verbose':
                    'Left Junctions<br>(Reference &rarr; Contig)'},
            {'field': 'right_junctions_html', 'verbose':
                    'Right Junctions<br>(Reference &rarr; Contig)'}
        ]


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

    # User specified data fields corresponding to the sample.
    # Examples: Growth rate, GFP amount, phenotype, # of mage cycles, etc.
    data = PostgresJsonField()

    # parent/child relations to other samples
    children = models.ManyToManyField('self',
        through='ExperimentSampleRelation',
        symmetrical=False,
        related_name='parents')

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

    def add_child(self, sample):
        """
        Create a relationship with another sample as as child.

        TODO: For now, don't complain if this is a parent sample as well,
        since we aren't doing anything fancy.

        Return True if successful.
        """
        return ExperimentSampleRelation.objects.get_or_create(
            parent= self,
            child= sample)

    def remove_child(self, sample):
        """
        Remove a parent/child relationship with another sample.

        Return True if present and removed.
        """
        child_relation = ExperimentSampleRelation.objects.filter(
            parent=self,
            child=sample)

        if child_relation:
            child_relation.delete()
            return True
        else:
            return False

    def get_children(self):
        """
        Use relationship table to get all children.
        """
        return self.children.all()

    def get_parents(self):
        """
        Use relationship table to get all parents.
        """
        return self.parents.all()

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

    @property
    def fastqc_links(self):
        """
        Links to the FASTQC output files.
        First checks if datasets are present, skips if missing.
        """

        links = []

        fqc_dataset_types = enumerate([
                Dataset.TYPE.FASTQC1_HTML,
                Dataset.TYPE.FASTQC2_HTML], start=1)

        for read_num, fqc_dataset_type in fqc_dataset_types:

            fastqc_dataset = get_dataset_with_type(self, fqc_dataset_type)

            if not fastqc_dataset:
                continue

            links.append(
                    '<a href="{url}" target="_blank">'
                    'Read {read_num}</a>'.format(
                            url=reverse(
                                    'main.views.fastqc_view',
                                    args=(self.project.uid, self.uid,
                                            read_num)),
                            read_num=read_num))

        return ', '.join(links)


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
            {'field': 'uid', 'verbose': 'Internal ID'},
            {'field': 'fastqc_links', 'verbose': 'FastQC'},
        ]


class ExperimentSampleRelation(UniqueUidModelMixin):
    """
    Explicit table linking parent and child samples.
    """
    parent = models.ForeignKey(ExperimentSample, related_name='parent_relationships')
    child = models.ForeignKey(ExperimentSample, related_name='child_relationships')


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

    def default_alignment_options():
        """
        Return the default alignment options.

        Includes currently:
        call_as_haploid : haploid calling mode (defaults to diploid)
        skip_het_only : remove het-only calls in diploid mode (default false)

        To do at some point:
        * custom arguments to bwa, gatk, freebayes, etc
        * enabling/changing of proecssing steps (DEFAULT_PROCESSING_MASK)
        """
        return json.dumps({
            'call_as_haploid': False,
            'skip_het_only': False
        })

    # see default_alignment_options()
    alignment_options = PostgresJsonField(default=default_alignment_options)

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

    # Statuses that indicate the alignment pipeline is running.
    PIPELINE_IS_RUNNING_STATUSES = [
        STATUS.ALIGNING,
        STATUS.VARIANT_CALLING
    ]

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

    def delete_model_data_dir(self):
        """Removes all data associated with this model.

        WARNING: Be careful with this method!
        """
        data_dir = self.get_model_data_dir()
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)

    @property
    def href(self):
        """Link to url view for this model.
        """
        return reverse(
                'main.views.alignment_view',
                args=(self.reference_genome.project.uid, self.uid))

    @property
    def run_time(self):
        """Time elapsed since alignment start.

        NOTE: This might be complicated by the not-so-clean implementation of
        the pipeline runner.
        """
        # Cases where alignment has not been run before.
        if (self.start_time is None or
                self.status == AlignmentGroup.STATUS.NOT_STARTED or
                self.status == AlignmentGroup.STATUS.UNKNOWN):
            return 'NOT RUNNING'

        # Determine effective end time to use for calculating running time,
        # depending on whether pipeline completed or not.
        if self.end_time is None:
            # Start time but no end time which typically should mean that
            # the pipeline is still running.

            # However, we still check for weird states because the pipeline
            # occasionally has issues.
            if self.status in [
                    AlignmentGroup.STATUS.FAILED,
                    AlignmentGroup.STATUS.COMPLETED]:
                return 'ERROR'

            effective_end_time = datetime.now()
        else:
            # End time exists so pipeline ran to completion or controlled
            # failure.
            effective_end_time = self.end_time

        # Return time delta, properly formatted.
        return re.match('(.*:.*:.*)\.',
                str(effective_end_time - self.start_time)).group(1)

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
                {'field':'end_time'},
                {'field':'run_time'}]

    def get_samples(self):
        """Many different tasks require getting the sample (or their UIDs)
        that are in this alignment group.
        """
        return ExperimentSample.objects.filter(
                experimentsampletoalignment__alignment_group=self)

    def get_or_create_vcf_output_dir(self):
        """Returns path to vcf root dir.
        """
        vcf_dir = os.path.join(self.get_model_data_dir(), 'vcf')
        ensure_exists_0775_dir(vcf_dir)
        return vcf_dir

    def get_combined_error_log_data(self):
        """Returns raw string representing entire error log for alignment.
        """
        vcf_dir = self.get_or_create_vcf_output_dir()

        # TODO(gleb): Support other error files.
        error_file = os.path.join(vcf_dir, 'merge_variant_data.error')

        if os.path.exists(error_file):
            with open(error_file) as fh:
                raw_data = fh.read()
        else:
            raw_data = 'None'
        return raw_data


class ExperimentSampleToAlignment(UniqueUidModelMixin):
    """Model that describes the alignment of a single ExperimentSample
    to an AlignmentGroup.
    """

    alignment_group = models.ForeignKey('AlignmentGroup')

    experiment_sample = models.ForeignKey('ExperimentSample')

    dataset_set = models.ManyToManyField('Dataset', blank=True, null=True)

    data = PostgresJsonField()

    class ASSEMBLY_STATUS:
        """
        The status of an Assembly
        """
        QUEUED = 'QUEUED TO ASSEMBLE'
        ASSEMBLING = 'ASSEMBLING'
        COMPLETED = 'COMPLETED'
        FAILED = 'FAILED'

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

    chromosome = models.ForeignKey('Chromosome')

    position = models.BigIntegerField('Position')

    ref_value = models.TextField('Ref')

    # User specified data fields corresponding to the variant
    data = PostgresJsonField()

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
    def variant_specific_tracks(self):
        return self.data.get(
                'variant_specific_tracks',
                {'alignment': [], 'coverage': []})

    @property
    def jbrowse_link(self):
        ref_genome_jbrowse = self.reference_genome.get_client_jbrowse_link()
        location_param = '&loc=' + str(self.position)
        full_href = ref_genome_jbrowse + location_param
        return '<a href="' + full_href + '">jbrowse</a>'

    @classmethod
    def get_field_order(clazz, **kwargs):
        raise NotImplementedError(
                "Currently, Variants are displayed via model_views.py")


class VariantCallerCommonData(Model, VisibleFieldMixin):
    """Model that describes data provided by a specific caller about a
    particular Variant.

    The reason for this model is that most variant callers are run for multiple
    ExperientSamples at the same time, generating some common data for each
    variant found, as well as data unique to each ExperimentSample. This model
    represents the common shared data.

    To be even more specific, the VCF format typically gives a row for each
    variant, where the first several columns describe the variant in general.
    This common data is stored in this model. There are additional columns in
    the vcf, one per ExperimentSample, which provides data about the
    relationship between the Variant and the ExperimentSample for that column.
    This data is stored in VariantEvidence instances, one per column.
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
        alt_value = self.alt_value
        if len(self.alt_value) > 10:
            alt_value = alt_value[:10] + '...'
        return 'var: ' + str(self.variant) + ', alt:' + alt_value

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


class SavedVariantFilterQuery(UniqueUidModelMixin):
    """Saved query belonging to the user.
    """
    owner = models.ForeignKey('UserProfile')

    text = models.TextField()


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


def get_or_create_derived_bam_dataset(sample_alignment, dataset_type,
        derivation_fn, force_rerun=False):
    """Gets or creates a new bam Dataset derived according to a provided function.

    The purpose of this function is to abstract the boilerplate that goes into
    creating a derived bam Dataset.

    Args:
        sample_alignment: ExperimentSampleToAlignment that is in a READY state.
        dataset_type: Dataset.TYPE of the dataset to get.
        derivation_fn: Function(sample_alignment, new_dataset).
            Mutates new_dataset. Should raise CalledProcessError if there is a
            problem during computing

    Returns:
        New Dataset.
    """
    # First, ensure the Dataset exists.
    new_dataset = get_dataset_with_type(
            sample_alignment, dataset_type)
    if new_dataset is None:
        new_dataset = Dataset.objects.create(
            label=dataset_type,
            type=dataset_type,
            status=Dataset.STATUS.NOT_STARTED)
        sample_alignment.dataset_set.add(new_dataset)

    # Next, check if the Dataset is already computed and can just be returned.
    if (not force_rerun and new_dataset.status == Dataset.STATUS.READY and
            os.path.exists(new_dataset.get_absolute_location())):
        return new_dataset

    # If here, we are going to run or re-run the Dataset so we reset the status
    # to indicate incomplete state.
    new_dataset.status = Dataset.STATUS.NOT_STARTED
    new_dataset.save(update_fields=['status'])

    try:
        # Start computing.
        new_dataset.status = Dataset.STATUS.COMPUTING
        new_dataset.save(update_fields=['status'])

        derivation_fn(sample_alignment, new_dataset)

        # Mark success.
        new_dataset.status = Dataset.STATUS.READY

    except subprocess.CalledProcessError:
        new_dataset.filesystem_location = ''
        new_dataset.status = Dataset.STATUS.FAILED

    new_dataset.save()

    return new_dataset
