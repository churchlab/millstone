"""
Signal registration.

See: https://docs.djangoproject.com/en/dev/topics/signals/.
"""

import pickle

from django.db.models.signals import m2m_changed
from django.db.models.signals import post_save

from models import AlignmentGroup
from models import Dataset
from models import ExperimentSample
from models import ExperimentSampleToAlignment
from models import ReferenceGenome
from models import VariantAlternate
from models import VariantEvidence
from models import VariantSet
from scripts.alignment_pipeline import ensure_bwa_index
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
from scripts.import_util import generate_fasta_from_genbank
from scripts.import_util import get_dataset_with_type
from scripts.jbrowse_util import prepare_jbrowse_ref_sequence
from scripts.snpeff_util import build_snpeff


# When a new ReferenceGenome is created, create its data dir.
def post_ref_genome_create(sender, instance, created, **kwargs):
    """Upon creation, create necessary data for jbrowse and snpeff."""
    if created:
        instance.ensure_model_data_dir_exists()
        instance.ensure_snpeff_dir()
        instance.ensure_jbrowse_dir()
# Run post-save commands after making a new ref genome object
post_save.connect(post_ref_genome_create, sender=ReferenceGenome,
        dispatch_uid='ref_genome_create')

def post_add_seq_to_ref_genome(sender, instance, **kwargs):
    """ When a dataset gets added to a ref_genome, we need to do some things
    if the ref_genome is a fasta or a genbank - run snpeff, jbrowse prep, and
    convert genbank to fasta. Each pk in pk_set is a Dataset object."""

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

    #Initialize variant key map field
    initialize_filter_key_map(instance)
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
    if not created or not 'gt_bases' in instance.data:
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


def post_sample_align_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_sample_create, sender=ExperimentSampleToAlignment,
        dispatch_uid='post_sample_create')


# We'll store freebayes and such under this location.
def post_alignment_group_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()
post_save.connect(post_alignment_group_create, sender=AlignmentGroup,
        dispatch_uid='alignment_group_create')
