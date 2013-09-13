"""
Methods for working with snpEff.
"""

from django import template
from Bio import SeqIO
import os
import os.path
import subprocess
import sys

from main.models import clean_filesystem_location
from main.models import Project
from main.models import ReferenceGenome
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.models import get_dataset_with_type
from scripts.import_util import sanitize_record_id
import settings

# TODO: These should be set somewhere else. snpeff_util and vcf_parser also use
#   them, but where should they go? settings.py seems logical, but it cannot
#   import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF


def build_snpeff(ref_genome, **kwargs):
    '''
    SnpEFF needs a config file for every reference genome, which lists a
    single reference genome, its chromosomes, and the codon table that
    each uses. For now we can assume that all our genomes will use bacterial
    codons. Every reference genome in the config file should look similar to:

    # Ecoli K12 MG1655
    NC_000913.genome : Escherichia_coli
        NC_000913.chromosomes : NC_000913
        NC_000913.NC_000913.codonTable: Bacterial_and_Plant_Plastid

    We have made a template that can do this with yaml rendering, in the
    snpEFF tools directory. Given a ref_genome object, it generates a
    snpEFF config file and builds and snpEFF database file for the genome,
    and places it in the ref genome's data dir under ./snpeff.
    '''

    # if no genbank file for this ref genome, then do nothing

    if not ref_genome.is_annotated():
        print "Snpeff indexing failed: No genbank for reference genome %s" % (
                ref_genome.uid)
        return

    #path to the snpeff dir for this project
    ref_genome.ensure_snpeff_dir()

    # Path we give to snp_eff should be the ref_genomes/uid/snpeff dir
    # and ref_genomes/uid/snpeff/uid is where snpeff will look for genes.gbk
    snpeff_uid_path = ref_genome.get_snpeff_directory_path()
    snpeff_path = os.path.join(get_snpeff_config_path(ref_genome))

    # Template data we will pass to the django template parser in order to
    # build the config file.
    templ_data = {}

    templ_data['snpeff_dir'] = snpeff_path

    # Put a soft link to the reference genome genbank file in the snpeff dir
    # under the ./snpeff/uid dir called genes.gb. Remove any old links
    # if present.
    ref_genome_path = ref_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()
    assert ref_genome_path is not None, "No reference source genbank."

    snpeff_genbank_symlink = os.path.join(snpeff_uid_path,'genes.gb')

    # Unlink if there was a link and then create a new link.
    try:
        os.unlink(snpeff_genbank_symlink)
    except OSError:
        # There was no symlink. That's fine.
        pass
    # Re-create the link.
    os.symlink(ref_genome_path, snpeff_genbank_symlink)

    # Fill in uid and chromosome data
    templ_data['uid'] = ref_genome.uid
    templ_data['chromosomes'] = []
    templ_data['label'] = ref_genome.label


    # Each record is a chromosome in the ref genome
    for seq_record in SeqIO.parse(
            open(ref_genome_path, "r"), "genbank"):

        # Add this record as a chromosome to this ref genome
        # TODO: Do we want to check seqrecords for sane/sanitized names?
        templ_data['chromosomes'].append(seq_record.id)

    # Render snpEff config template
    render_snpeff_config(templ_data, os.path.join(snpeff_path,'snpeff.config'))

    # Build snpEff database
    build_snpeff_db(os.path.join(snpeff_path,'snpeff.config'), ref_genome.uid)


def render_snpeff_config(
        data,
        output_file_location,
        config_template_location= settings.SNPEFF_CFG_TEMPLATE_PATH):
    """Renders the config file given the data.
    Args:
        data: Dictionary of data where keys are the template variables
            and values are the data to render into the template.
    """
    with open(config_template_location) as template_fh:
        config_template = template.Template(template_fh.read())

    context = template.Context(data)

    with open(output_file_location, 'w') as output_fh:
        output_fh.write(config_template.render(context))

def build_snpeff_db(snpeff_config_path, ref_genome_uid):
    """
    Call snpeff's build function to build a database based on the
    genbank file. Command should look like:

    java -jar snpEff.jar build -genbank -v $UID -c $CFG_FILE
    """

    snpeff_args = [
        'java',
        '-jar', settings.SNPEFF_JAR_PATH,
        'build',
        '-genbank',
        '-v', ref_genome_uid,
        '-c', snpeff_config_path,
        '-q',
        '-noLog'
    ]

    # TODO: this redirect breaks nose tests, so
    # # If we need to debug the build step, don't throw away the output here
    # if settings.SNPEFF_BUILD_DEBUG:
    #     stdout_dest = sys.__stdout__
    #     stderr_dest = sys.__stderr__
    # else:
    #     stdout_dest=open(os.devnull, 'w')
    #     stderr_dest=open(os.devnull, 'w')

    snpeff_proc = subprocess.call(snpeff_args)

if __name__ == '__main__':
    gbk_test_path = ['/scratch/dbg/bw_strains/gbk/BW25113_updated.gbk']
    render_snpeff_config(
        generate_snpeff_config(
            gbk_test_path,
            ['BW25113_ref_genome'],
            '/path/to/project'),
        '/scratch/dbg/bw_strains/gbk/BW25113_updated.snpeff.config')

def run_snpeff(alignment_group, alignment_type):
    """Run snpeff on an alignment group after creating a vcf with a snpcaller.

    We only use the alignment type to store the snpeff file.
    """

    # Get the reference genome uid to get the config path and snpeff genome name
    ref_genome = alignment_group.reference_genome
    ref_genome_uid = alignment_group.reference_genome.uid
    snpeff_config_path = get_snpeff_config_path(ref_genome)

    # Get the freebayes vcf file as input
    assert alignment_group.dataset_set.get(type=VCF_DATASET_TYPE)
    vcf_input_filename = alignment_group.dataset_set.get(
            type=VCF_DATASET_TYPE).get_absolute_location()

    # Prepare a directory to put the output file.
    # We'll put them in /projects/<project_uid>/alignment_groups/vcf/snpeff/
    #     <alignment_type>.vcf
    # We'll save these for now, maybe it's not necessary later.

    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    snpeff_vcf_dir = os.path.join(vcf_dir, 'snpeff')
    ensure_exists_0775_dir(snpeff_vcf_dir)
    vcf_output_filename = os.path.join(
            snpeff_vcf_dir, alignment_type + '.vcf')

    snpeff_args = [
        'java',
        '-jar', settings.SNPEFF_JAR_PATH,
        'eff',
        '-v',
        '-i',
        'vcf',
        '-o',
        'vcf',
        '-c', os.path.join(get_snpeff_config_path(ref_genome),'snpeff.config'),
        '-ud', str(settings.SNPEFF_UD_INTERVAL_LENGTH),
#        '-q',
#        '-noLog',
#        '-t', str(settings.SNPEFF_THREADS),
        ref_genome_uid,
        vcf_input_filename
    ]

    print snpeff_args

    with open(vcf_output_filename, 'w') as fh_out:
        snpeff_proc = subprocess.call(
            snpeff_args,
            stdout=fh_out)

    #Add the snpeff VCF file to the dataset_set
    # If a Dataset already exists, delete it, might have been a bad run.
    existing_set = Dataset.objects.filter(
            type=VCF_ANNOTATED_DATASET_TYPE,
            label=VCF_ANNOTATED_DATASET_TYPE,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    if len(existing_set) > 0:
        existing_set[0].delete()

    dataset = Dataset.objects.create(
            type=VCF_ANNOTATED_DATASET_TYPE,
            label=VCF_ANNOTATED_DATASET_TYPE,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    alignment_group.dataset_set.add(dataset)

def get_snpeff_config_path(ref_genome):
    """The parent directory of the snpeff model dir holds the config file, but
    the model creates another directory under that that actually holds all of
    the reference genome's snpeff data.

    This function returns that parent dir that holds the per-ReferenceGenome
    config file.
    """
    return os.path.abspath(os.path.join(ref_genome.get_snpeff_directory_path(),
            os.pardir))