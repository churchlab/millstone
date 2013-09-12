from django.db.models.signals import post_save
from django.db.models.signals import m2m_changed
from scripts.jbrowse_util import prepare_jbrowse_ref_sequence
from scripts.snpeff_util import build_snpeff
from scripts.import_util import generate_fasta_from_genbank

from models import ReferenceGenome
from models import Dataset

# When a new ReferenceGenome is created, create its data dir.
def post_ref_genome_create(sender, instance, created, **kwargs):
    """ Upon creation, create necessary data for jbrowse and snpeff."""
    if created:
        instance.ensure_model_data_dir_exists()
        instance.ensure_snpeff_dir()
        instance.ensure_jbrowse_dir()

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

# Run post-save commands after making a new ref genome object
post_save.connect(post_ref_genome_create, sender=ReferenceGenome,
        dispatch_uid='ref_genome_create')

# Run snpeff and jbrowse housekeeping after linking seq file dataset to a
#   reference genome obj
m2m_changed.connect(post_add_seq_to_ref_genome,
        sender=ReferenceGenome.dataset_set.through,
        dispatch_uid='add_seq_to_ref_genome')

