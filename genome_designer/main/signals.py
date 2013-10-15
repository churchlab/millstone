from django.db.models.signals import post_save
from django.db.models.signals import m2m_changed
from scripts.jbrowse_util import prepare_jbrowse_ref_sequence
from scripts.snpeff_util import build_snpeff
from scripts.import_util import generate_fasta_from_genbank
from scripts.dynamic_snp_filter_key_map import initialize_filter_key_map
import pickle

from models import ReferenceGenome
from models import Dataset
from models import VariantSet
from models import VariantAlternate
from models import VariantEvidence

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

    #Initialize variant key map field
    initialize_filter_key_map(instance)

def post_variant_evidence_create(sender, instance, created, **kwargs):
    """Add existing VariantAlternates to this VariantEvidence Object."""

    if not created or not 'gt_bases' in instance.data:
        return

    gt_bases = pickle.loads(instance.data['gt_bases'])

    #if this variant evidence is a non-call, no need to add alternate alleles
    if gt_bases is None: return

    assert ('|' not in gt_bases), (
            'GT bases string is phased;' + 
            'this is not handled and should never happen...')

    gt_bases = gt_bases.split('/')

    for gt_base in gt_bases:
        
        try:
            variant = instance.variant_caller_common_data.variant

            #Skip if this is not an alternate allele            
            if variant.ref_value == gt_base: continue

            instance.variantalternate_set.add(VariantAlternate.objects.get(
                    variant=variant,
                    alt_value=gt_base
                ))


        except VariantAlternate.DoesNotExist:
            print ('Attempt to add a SampleEvidence with an alternate ' +
                    'allele that is not present for this variant!')
            raise

def post_variant_set_create(sender, instance, created, **kwargs):
    if created:
        instance.ensure_model_data_dir_exists()

# Run post-save commands after making a new ref genome object
post_save.connect(post_ref_genome_create, sender=ReferenceGenome,
        dispatch_uid='ref_genome_create')

# Run post-save commands after making a new variant evidence object
post_save.connect(post_variant_evidence_create, sender=VariantEvidence,
        dispatch_uid='variant_evidence_create')

# Run snpeff and jbrowse housekeeping after linking seq file dataset to a
#   reference genome obj
m2m_changed.connect(post_add_seq_to_ref_genome,
        sender=ReferenceGenome.dataset_set.through,
        dispatch_uid='add_seq_to_ref_genome')

# When a new VariantSet is created, create the data directory.
post_save.connect(post_variant_set_create, sender=VariantSet,
        dispatch_uid='variant_set_create')


