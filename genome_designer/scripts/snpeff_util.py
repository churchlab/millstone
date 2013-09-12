"""
Methods for working with snpEff.
"""

from django import template
import settings
from Bio import SeqIO
import os
import os.path
import subprocess

from main.models import Project
from main.models import ReferenceGenome
from main.models import Dataset


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

    if not ref_genome.dataset_set.filter(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).exists():
        print "Snpeff indexing failed: No genbank for reference genome %s" % (
                ref_genome.uid)
        return

    #path to the snpeff dir for this project
    ref_genome.ensure_snpeff_dir()

    # Path we give to snp_eff should be the ref_genomes/uid/snpeff dir
    # and ref_genomes/uid/snpeff/uid is where snpeff will look for genes.gbk
    snpeff_uid_path = ref_genome.get_snpeff_directory_path()
    snpeff_path = os.path.abspath(os.path.join(snpeff_uid_path, os.pardir))

    # Template data we will pass to the django template parser in order to
    # build the config file.
    templ_data = {}

    templ_data['snpeff_dir'] = snpeff_path

    #put a soft link to the reference genome genbank file in the snpeff dir
    # under the ./snpeff/uid dir called genes.gb
    ref_genome_path = ref_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()
    assert ref_genome_path is not None, "No reference soure genbank."
    os.symlink(ref_genome_path, os.path.join(snpeff_uid_path,'genes.gb'))

    # Fill in uid and chromosome data
    templ_data['uid'] = ref_genome.uid
    templ_data['chromosomes'] = []

    # Each record is a chromosome in the ref genome
    for seq_record in SeqIO.parse(
            open(ref_genome_path, "r"), "genbank"):

        # Add this record as a chromosome to this ref genome
        # TODO: Do we want to check seqrecords for sane/sanitized names?
        templ_data['chromosomes'].append(seq_record.name)

    # Render snpEff config template
    snpeff_config_file = os.path.join(snpeff_path, 'snpeff.config')
    render_snpeff_config(templ_data, snpeff_config_file)

    # Build snpEff database
    build_snpeff_db(snpeff_config_file, ref_genome.uid)


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
        '-jar',
        settings.SNPEFF_JAR_PATH,
        'build',
        '-genbank',
        '-v',
        ref_genome_uid,
        '-c',
        snpeff_config_path,
        '-q',
        '-noLog'
    ]
    snpeff_proc = subprocess.call(snpeff_args,
        stdout=open(os.devnull, 'w'),
        stderr=open(os.devnull, 'w'))

if __name__ == '__main__':
    gbk_test_path = ['/scratch/dbg/bw_strains/gbk/BW25113_updated.gbk']
    render_snpeff_config(
        generate_snpeff_config(
            gbk_test_path,
            ['BW25113_ref_genome'],
            '/path/to/project'),
        '/scratch/dbg/bw_strains/gbk/BW25113_updated.snpeff.config')


