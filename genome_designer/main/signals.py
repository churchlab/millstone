"""
Signal registration.

See: https://docs.djangoproject.com/en/dev/topics/signals/.
"""

import os
import shutil

from Bio import SeqIO
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User
from django.db.models.signals import m2m_changed
from django.db.models.signals import post_save
from django.db.models.signals import pre_delete
from django.db.models.signals import post_delete

from models import AlignmentGroup
from models import Contig
from models import Dataset
from models import ExperimentSample
from models import ExperimentSampleToAlignment
from models import Project
from models import ReferenceGenome
from models import UserProfile
from models import VariantEvidence
from models import VariantSet
from models import VariantToVariantSet
from utils.import_util import prepare_ref_genome_related_datasets
from utils.import_util import add_chromosomes
from utils.import_util import sanitize_sequence_dataset
from variants.dynamic_snp_filter_key_map import initialize_filter_key_map
from variants.dynamic_snp_filter_key_map import update_sample_filter_key_map


# Since the registration flow creates a django User object, we want to make
# sure that the corresponding UserProfile is also created
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(user=instance)
post_save.connect(create_user_profile, sender=User,
        dispatch_uid='user_profile_create')


# Delete the associated Django User on UserProfile delete.
def post_user_profile_delete(sender, instance, **kwargs):
    # This is only necessary if the UserProfile was deleted directly.
    # If the User model was deleted, then the User will no longer exist
    # by the time we get here, so nothing to do.
    user_profile = instance
    maybe_user_matches = User.objects.filter(id=user_profile.user_id)
    if len(maybe_user_matches):
        assert len(maybe_user_matches) == 1
        user = maybe_user_matches[0]
        assert user.user_profile_id == user_profile.id
        maybe_user_matches[0].delete()
post_delete.connect(post_user_profile_delete, sender=UserProfile,
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


def post_contig_delete(sender, instance, **kwargs):
    """Delete all data associated with Contig including:
        * contig data
        * associated variant caller data.
    """
    data_dir = instance.get_model_data_dir()
    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)

    # In case vccd has already been deleted
    try:
        vccd = instance.variant_caller_common_data
    except (ObjectDoesNotExist, AttributeError):
        return
    if vccd is None:
        return

    # If this is last VariantCallerCommonData for Variant, delete the whole
    # Variant.
    vccd.variant.reference_genome.invalidate_materialized_view()
    if vccd.variant.variantcallercommondata_set.count() == 1:
        vccd.variant.delete()
    else:
        vccd.delete()
post_delete.connect(post_contig_delete, sender=Contig,
        dispatch_uid='contig_delete')


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
    instance.drop_materialized_view()
pre_delete.connect(pre_ref_genome_delete, sender=ReferenceGenome,
        dispatch_uid='ref_genome_delete')


def post_add_seq_to_ref_genome(sender, instance, **kwargs):
    """When a dataset gets added to a ref_genome, we need generate some
    additional data via prepare_ref_genome_related_datasets().

    If the ref_genome is a fasta or a genbank - run snpeff, jbrowse prep, and
    convert genbank to fasta and gff. Each pk in pk_set is a Dataset object.

    NOTE: If the associated Dataset has the NOT_READY status, then this signal
    will abort before doing anything interesting and it's up to the client
    to handle calling prepare_ref_genome_related_datasets().
    """
    # Skip unless we've already added the relationship
    if kwargs['action'] != 'post_add':
        return

    # TODO: Why is this a loop? When will this ever have more than one pk in
    # kwargs['pk_set']?
    model = kwargs['model']
    for pk in kwargs['pk_set']:
        dataset = model.objects.get(pk=pk)
        if dataset.status == Dataset.STATUS.NOT_STARTED:
            continue
        sanitize_sequence_dataset(dataset)
        prepare_ref_genome_related_datasets(instance, dataset)
        add_chromosomes(instance, dataset)

m2m_changed.connect(post_add_seq_to_ref_genome,
    sender=ReferenceGenome.dataset_set.through,
    dispatch_uid='add_seq_to_ref_genome')


def post_add_seq_to_contig(sender, instance, **kwargs):
    """When a dataset gets added to a contig, we need to grab
    the number of bases from the dataset if it is a fasta and
    ensure that it is not a multifasta
    """

    # Skip unless we've already added the relationship
    if kwargs['action'] != 'post_add':
        return

    # TODO: Why is this a loop? When will this ever have more than one pk in
    # kwargs['pk_set']?
    model = kwargs['model']
    for pk in kwargs['pk_set']:
        dataset = model.objects.get(pk=pk)
        if dataset.type == Dataset.TYPE.REFERENCE_GENOME_FASTA:
            assert instance.num_bases == 0
            seq_record_count = 0
            for seq_record in SeqIO.parse(dataset.get_absolute_location(),
                    'fasta'):
                if seq_record_count != 0:
                    raise Exception('Contig associated fastas can only' +
                                    'contain one sequence')
                instance.num_bases = len(seq_record)
                seq_record_count += 1

m2m_changed.connect(post_add_seq_to_contig,
        sender=Contig.dataset_set.through,
        dispatch_uid='post_add_seq_to_contig')


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


def pre_alignment_group_delete(sender, instance, **kwargs):
    instance.delete_model_data_dir()
    instance.reference_genome.drop_materialized_view()
pre_delete.connect(pre_alignment_group_delete, sender=AlignmentGroup,
        dispatch_uid='alignment_group_delete')


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
