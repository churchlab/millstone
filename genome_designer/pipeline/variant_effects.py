"""
Methods for working with snpEff.
"""

from django import template
from Bio import SeqIO
import os
import os.path
import subprocess
import sys
import vcf
from collections import defaultdict
from collections import OrderedDict
import re
from itertools import chain
from string import Template
from StringIO import StringIO

from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.models import get_dataset_with_type
from main.model_utils import clean_filesystem_location
from utils import ensure_line_lengths
from utils import uppercase_underscore
import settings

# TODO: These should be set somewhere else. snpeff_util and vcf_parser also use
#   them, but where should they go? settings.py seems logical, but it cannot
#   import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF

#Compile this SNPEFF parsing refex only once
# Example:
#   NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC||CODING
#        |b3038|1|1)

SNPEFF_FIELDS = OrderedDict([
    ('EFFECT', {
        'id':'EFFECT',
        'type':'String',
        'description':'Effect type of Variant.'}
    ),
    ('IMPACT', {
        'id':'IMPACT',
        'type':'String',
        'description':'Effect impact {High, Moderate, Low, Modifier}.'}
    ),
    ('CLASS', {
        'id':'CLASS',
        'type':'String',
        'description':'Functional class {NONE, SILENT, MISSENSE, NONSENSE}.'}
    ),
    ('CONTEXT', {
        'id':'CONTEXT',
        'type':'String',
        'description':'old_codon/new_codon OR distance to transcript.'}
    ),
    ('AA', {
        'id':'AA',
        'type':'String',
        'description':'Amino acid change: old_AA AA_position/new_AA.'}
    ),
    ('TRLEN', {
        'id':'TRLEN',
        'type':'Integer',
        'description':'Length of protein in amino acids.'}
    ),
    ('GENE', {
        'id':'GENE',
        'type':'String',
        'description':'Gene Name.'}
    ),
    ('BIOTYPE', {
        'id':'BIOTYPE',
        'type':'String',
        'description':'Transcript bioType, if available.'}
    ),
    ('CODING', {
        'id':'CODING',
        'type':'String',
        'description':'Either CODING or NONCODING.'}
    ),
    ('TR', {
        'id':'TR',
        'type':'String',
        'description':'Transcript ID (usually ENSEMBL IDs).'}
    ),
    ('RANK', {
        'id':'RANK',
        'type':'String',
        'description':'Exon rank or Intron rank.'}
    ),
    ('GT', {
        'id':'GT',
        'type':'String',
        'description':'Genotype number corresponding to this effect.'}
    ),
    ('ERR', {
        'id':'ERR',
        'type':'String',
        'description':'Any Errors.'}
    ),
    ('WARN', {
        'id':'WARN',
        'type':'String',
        'description':'Any Warnings.'}
    )]
)

SNPEFF_INFO_TEMPLATE = Template(','.join([
        '##INFO=<ID=EFF_$id',
        'Number=A',
        'Type=$type',
        'Description="$description">']))

SNPEFF_ALT_RE = re.compile(r''.join([
        r'(?P<{:s}>\w+)\((?P<{:s}>[^\|]*)',
        r'\|(?P<{:s}>[^\|]*)' * (len(SNPEFF_FIELDS.keys())-4),
        r'\|?(?P<{:s}>[^\|]*)\|?(?P<{:s}>[^\|]*)\)']
        ).format(*SNPEFF_FIELDS.keys()))

def build_snpeff(ref_genome):
    """
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
    """

    # if no genbank file for this ref genome, then do nothing
    if not ref_genome.is_annotated():
        print "Snpeff indexing failed: No genbank for reference genome %s" % (
                ref_genome.uid)
        return

    # Path we give to snp_eff should be the ref_genomes/uid/snpeff dir
    # and ref_genomes/uid/snpeff/uid is where snpeff will look for genes.gbk
    ref_genome.ensure_snpeff_dir()
    snpeff_uid_path = ref_genome.get_snpeff_directory_path()
    snpeff_path = os.path.join(get_snpeff_config_path(ref_genome))

    # Template data we will pass to the django template parser in order to
    # build the config file.
    templ_data = {}

    templ_data['snpeff_dir'] = snpeff_path

    # Put a soft link to the reference genome genbank file in the snpeff dir
    # under the ./snpeff/uid dir called genes.gb. Remove any old links
    # if present.
    ref_genome_path = get_dataset_with_type(ref_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()
    assert ref_genome_path is not None, "No reference source genbank."

    snpeff_genbank_filename = os.path.join(snpeff_uid_path,'genes.gb')

    # TODO: The symlink below doesn't work - we need to now directly modify
    # the Genbank for it to work with SnpEFF, so we'll keep a modified copy in
    # the Ref Genome's SnpEFF directory. We have to modify it to ensure that
    # the name and ID are the same and to ensure minimum line lengths.

#    # Unlink if there was a link and then create a new link.
#    try:
#        os.unlink(snpeff_genbank_filename)
#    except OSError:
#        # There was no symlink. That's fine.
#        pass
#    # Re-create the link.
#    os.symlink(ref_genome_path, snpeff_genbank_filename)

    # Fill in uid and chromosome data
    templ_data['uid'] = ref_genome.uid
    templ_data['chromosomes'] = []
    templ_data['label'] = ref_genome.label

    new_gb = []

    # Each record is a chromosome in the ref genome
    for seq_record in SeqIO.parse(
            open(ref_genome_path, "r"), "genbank"):

        # Set the ACCESSION/LOCUS/VERSION to all be the same for this
        # new modified genbank
        seq_record.id = seq_record.name
        new_gb.append(seq_record)

        # Add this record as a chromosome to this ref genome
        # TODO: Do we want to check seqrecords for sane/sanitized names?
        templ_data['chromosomes'].append(seq_record.name)

    # Save a modified copy of the genbank for snpEff
    # with open(snpeff_genbank_symlink, 'w') as snpeff_gbk_file:
    SeqIO.write(new_gb, snpeff_genbank_filename, "genbank")
    # Stop-gap fix to ensure line lengths in GENbANK to appease SNPEFF
    ensure_line_lengths(snpeff_genbank_filename)

    templ_data['chrs_string'] = ','.join(templ_data['chromosomes'])

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

    print >> sys.stderr, ' '.join(snpeff_args)

    # TODO: this redirect breaks nose tests, so
    # # If we need to debug the build step, don't throw away the output here
    # if settings.SNPEFF_BUILD_DEBUG:
    #     stdout_dest = sys.__stdout__
    #     stderr_dest = sys.__stderr__
    # else:
    #     stdout_dest=open(os.devnull, 'w')
    #     stderr_dest=open(os.devnull, 'w')

    snpeff_proc = subprocess.call(snpeff_args)


def get_snpeff_vcf_output_path(alignment_group, alignment_type):
    """Returns the path to SnpEff dir, reltaive to the AlignmentGroup data
    location.

    Ensures that the intermediate dirs are created so the file can be written.

    Path is of the form:
    /projects/<project_uid>/alignment_groups/vcf/snpeff/<alignment_type>.vcf
    """
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    snpeff_vcf_dir = os.path.join(vcf_dir, 'snpeff')
    ensure_exists_0775_dir(snpeff_vcf_dir)
    vcf_output_filename = os.path.join(
            snpeff_vcf_dir, uppercase_underscore(alignment_type) + '.vcf')
    return vcf_output_filename


def run_snpeff(alignment_group, alignment_type):
    """Run snpeff on an alignment group after creating a vcf with a snpcaller.

    We only use the alignment type to store the snpeff file.

    Returns the snpeff vcf output filename.
    """

    # Get the reference genome uid to get the config path and snpeff genome name
    ref_genome = alignment_group.reference_genome
    ref_genome_uid = alignment_group.reference_genome.uid
    snpeff_config_path = get_snpeff_config_path(ref_genome)

    # Get the freebayes vcf file as input
    assert get_dataset_with_type(alignment_group,
            type=VCF_DATASET_TYPE)
    vcf_input_filename = get_dataset_with_type(alignment_group,
            type=VCF_DATASET_TYPE).get_absolute_location()

    # Prepare a directory to put the output file.
    vcf_output_filename = get_snpeff_vcf_output_path(alignment_group,
            alignment_type)

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
        '-q',
        '-noLog',
#        '-t', str(settings.SNPEFF_THREADS),
        ref_genome_uid,
        vcf_input_filename
    ]

    print ' '.join(snpeff_args)

    with open(vcf_output_filename, 'w') as fh_out:

        #Create a stringIO buffer to hold the unedited snpeff vcf
        raw_snpeff_vcf = StringIO()

        snpeff_proc = subprocess.Popen(
            snpeff_args,
            stdout=subprocess.PIPE)
        convert_snpeff_info_fields(snpeff_proc.stdout, fh_out)

    return vcf_output_filename

    #clean up the snpEff_genes.txt and snpEff_summary.txt files in the home dir
    for file in settings.SNPEFF_SUMMARY_FILES:
        os.remove(os.path.join(os.getcwd(),file))

def convert_snpeff_info_fields(vcf_input_fh, vcf_output_fh):
    """This function takes a VCF file on an input stream, reads it in,
    converts the single EFF field to a set of EFF fields, and then returns
    the modified VCF file on an output stream.

    The snpeff field starts out as a long string, consisting of many fields
    each separated by pipes.

    Effects information is added to the INFO field using an 'EFF' tag.
    There can be multiple effects separated by comma. The format for each
    effect is:

    Effect ( Effect_Impact | Codon_Change | Amino_Acid_change | Gene_Name
            | Gene_BioType | Coding | Transcript | Rank [ | ERRORS
            | WARNINGS ] )

    Details for each field are here:
        http://snpeff.sourceforge.net/SnpEff_manual.html

    We will pull out all of these fields separately into INFO_EFF_* and return
    a new VCF file.
    """

    vcf_reader = vcf.Reader(vcf_input_fh)

    # Generate extra header rows.
    # TODO: This method is internal to pyVCF, so if they change it,
    # this will break. Maybe we should copy their code?
    parser = vcf.parser._vcf_metadata_parser()
    for field, values in SNPEFF_FIELDS.items():

        # Create a new header line from the new field.
        print values
        new_header_line = SNPEFF_INFO_TEMPLATE.substitute(values)

        # Add this extra header line to the vcf reader.
        vcf_reader._header_lines.append(new_header_line)

        # Parse the header line as an Info obj, add it to the reader.
        key, val = parser.read_info(new_header_line)
        vcf_reader.infos[key] = val

    vcf_writer = vcf.Writer(vcf_output_fh, vcf_reader)

    # Write the old records with the new EFF INFO fields
    for record in vcf_reader:
        vcf_writer.write_record(populate_record_eff(record))

def populate_record_eff(vcf_record):
    """This function takes a single VCF record and separates the EFF key
    from snpEFF into a variety of individual keys.

    The snpeff field starts out as a long string, consisting of many fields
    each separated by pipes.

    Effects information is added to the INFO field using an 'EFF' tag.
    There can be multiple effects separated by comma. The format for each
    effect is:

    Effect ( Effect_Impact | Codon_Change | Amino_Acid_change | Gene_Name
            | Gene_BioType | Coding | Transcript | Rank [ | ERRORS
            | WARNINGS ] )

    Example:
    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC||CODING
        |b3038|1|1)

    Details for each field are here:
        http://snpeff.sourceforge.net/SnpEff_manual.html

    NOTE: It is possible to have multiple values for each field if there are
    multiple alts:

    Example for REF=T;ALT=[C,G]:
    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC||CODING
        |b3038|1|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aGg|M239T
        |386|ygiC||CODING|b3038|1|1)

    In the code below, the above example would yield two re.groupdict() objs
    from eff_group_iterators. These would be chained together and the values
    would be 'zipped', so that the EFF_AA field would be a list of two values:
    ['aTg/aCg','aTg/aGg'].
    """

    # Get the Eff field for this record.
    eff_field_lists = defaultdict(list)
    assert hasattr(vcf_record,'INFO'), 'No INFO attr, not a vcf record'
    assert 'EFF' in vcf_record.INFO, 'VCF record has no EFF INFO field'
    value = vcf_record.INFO['EFF']

    # Find iter produces a separate set of groups for every comma-separated
    # alt EFF field set
    eff_group_iterator = (SNPEFF_ALT_RE.match(i) for i in value)

    # This chains together a list of all EFF fields from all alts
    eff_fields = list(chain.from_iterable(
            (eff.groupdict().items() for eff in eff_group_iterator)))

    for k, v in eff_fields:
        # mark empty fields as 'bad'
        if v == '': v = '.'
        eff_field_lists['EFF_'+k].append(v)

    vcf_record.INFO.update(eff_field_lists)

    return vcf_record

def get_snpeff_config_path(ref_genome):
    """The parent directory of the snpeff model dir holds the config file, but
    the model creates another directory under that that actually holds all of
    the reference genome's snpeff data.

    This function returns that parent dir that holds the per-ReferenceGenome
    config file.
    """
    return os.path.abspath(os.path.join(ref_genome.get_snpeff_directory_path(),
            os.pardir))
