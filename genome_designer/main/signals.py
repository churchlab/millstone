"""
Signal registration.

See: https://docs.djangoproject.com/en/dev/topics/signals/.
"""

from django.contrib.auth.models import User
from django.db.models.signals import m2m_changed
from django.db.models.signals import post_save
from django.db.models.signals import pre_delete

from models import AlignmentGroup
from models import Dataset
from models import ExperimentSample
from models import ExperimentSampleToAlignment
from models import Project
from models import ReferenceGenome
from models import UserProfile
from models import VariantEvidence
from models import VariantSet
from models import VariantToVariantSet
from pipeline.read_alignment_util import ensure_bwa_index
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.dynamic_snp_filter_key_map import update_sample_filter_key_map
from scripts.import_util import generate_fasta_from_genbank
from scripts.import_util import generate_gff_from_genbank
from scripts.import_util import get_dataset_with_type
from scripts.jbrowse_util import prepare_jbrowse_ref_sequence
from scripts.jbrowse_util import add_genbank_file_track
from pipeline.variant_effects import build_snpeff


# Since the registration flow creates a django User object, we want to make
# sure that the corresponding UserProfile is also created
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(user=instance)
post_save.connect(create_user_profile, sender=User,
        dispatch_uid='user_profile_create')


# When a new Project is created, create the data directory.
def post_project_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_project_create, sender=Project,
        dispatch_uid='project_create')


# Delete all Project data when it is deleted.
def pre_project_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
pre_delete.connect(pre_project_delete, sender=Project,
        dispatch_uid='project_delete')


# When a new ReferenceGenome is created, create its data dir.
def post_ref_genome_create(sender, instance, created, **kwargs):
    """Upon creation, create necessary data for jbrowse and snpeff."""
    if created:
        instance.ensure_model_data_dir_exists()
        instance.ensure_snpeff_dir()
        instance.ensure_jbrowse_dir()


        # Key Map should start out empty...
        assert not any(instance.variant_key_map)

        instance.variant_key_map = initialize_filter_key_map()
        instance.save()
post_save.connect(post_ref_genome_create, sender=ReferenceGenome,
        dispatch_uid='ref_genome_create')


# Delete all Project data when it is deleted.
def pre_ref_genome_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
pre_delete.connect(pre_ref_genome_delete, sender=ReferenceGenome,
        dispatch_uid='ref_genome_delete')


def post_add_seq_to_ref_genome(sender, instance, **kwargs):
    """ When a dataset gets added to a ref_genome, we need to do some things
    if the ref_genome is a fasta or a genbank - run snpeff, jbrowse prep, and
    convert genbank to fasta and gff. Each pk in pk_set is a Dataset object."""

    #Skip unless we've already added the relationship
    if kwargs['action'] != 'post_add':
        return

    model = kwargs['model']
    for pk in kwargs['pk_set']:
        dataset = model.objects.get(pk=pk)
        if dataset.type == Dataset.TYPE.REFERENCE_GENOME_FASTA:
            # Run jbrowse ref genome processing
            prepare_jbrowse_ref_sequence(instance)

        elif dataset.type == Dataset.TYPE.REFERENCE_GENOME_GENBANK:
            # Run snpeff build after creating ReferenceGenome obj
            build_snpeff(instance)
            # If no FASTA exists, then create it from the genbank
            # (this is done post-save so it should work if both files are
            # present)
            if not instance.dataset_set.filter(
                    type=Dataset.TYPE.REFERENCE_GENOME_FASTA).exists():
                generate_fasta_from_genbank(instance)
                generate_gff_from_genbank(instance)

            # Run jbrowse genbank genome processing for genes
            add_genbank_file_track(instance)

    # Run snpeff and jbrowse housekeeping after linking seq file dataset to a
    #   reference genome obj

    # Create the bwa index before perfoming the alignments in parallel.
    ref_genome_fasta = get_dataset_with_type(
            instance,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    ensure_bwa_index(ref_genome_fasta)

m2m_changed.connect(post_add_seq_to_ref_genome,
    sender=ReferenceGenome.dataset_set.through,
    dispatch_uid='add_seq_to_ref_genome')


def post_variant_evidence_create(sender, instance, created, **kwargs):
    """Add existing VariantAlternates to this VariantEvidence Object."""
    if not created or not 'GT_BASES' in instance.data:
        return
    instance.create_variant_alternate_association()
post_save.connect(post_variant_evidence_create, sender=VariantEvidence,
        dispatch_uid='variant_evidence_create')


def post_variant_set_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
# When a new VariantSet is created, create the data directory.
post_save.connect(post_variant_set_create, sender=VariantSet,
        dispatch_uid='variant_set_create')


def post_sample_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_sample_create, sender=ExperimentSample,
        dispatch_uid='post_sample_create')


def pre_sample_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
pre_delete.connect(pre_sample_delete, sender=ExperimentSample,
        dispatch_uid='post_sample_delete')


def post_sample_align_create(sender, instance, created, **kwargs):
    '''
    Make a model data dir and update variant filter keys with
    the sample metadata field names.
    '''
    if created:
        instance.ensure_model_data_dir_exists()
        ref_genome = instance.alignment_group.reference_genome
        # Update from the database
        ref_genome = ReferenceGenome.objects.get(id=ref_genome.id)

        # Sample metadata field names
        if not any(ref_genome.variant_key_map):
            raise ValueError('Variant Key Map was not properly initialized'
                    ' for %s: %s' % (ref_genome.uid,
                            ref_genome.variant_key_map))

        experiment_sample = instance.experiment_sample
        update_sample_filter_key_map(ref_genome, experiment_sample)
post_save.connect(post_sample_align_create, sender=ExperimentSampleToAlignment,
        dispatch_uid='post_sample_align_create')


def post_alignment_group_create(sender, instance, created, **kwargs):
    '''
    We'll store freebayes and such under this location.
    '''
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_alignment_group_create, sender=AlignmentGroup,
        dispatch_uid='alignment_group_create')


def post_vtvs_save(sender, instance, created, **kwargs):
    instance.variant.reference_genome.invalidate_materialized_view()
post_save.connect(post_vtvs_save, sender=VariantToVariantSet,
        dispatch_uid='vtvs_save')


def post_vtvs_delete(sender, instance, **kwargs):
    instance.variant.reference_genome.invalidate_materialized_view()
pre_delete.connect(post_vtvs_delete, sender=VariantToVariantSet,
        dispatch_uid='vtvs_delete')


def post_sample_alignment_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
pre_delete.connect(post_sample_alignment_delete,
        sender=ExperimentSampleToAlignment,
        dispatch_uid='sample_alignment_delete')


def post_variant_set_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_variant_set_create, sender=VariantSet,
        dispatch_uid='variant_set_create')
