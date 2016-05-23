"""
Methods for working with snpEff.
"""

from collections import defaultdict
from collections import OrderedDict
from itertools import chain
import os
import re
from string import Template
import subprocess
import sys

from Bio import SeqIO
from django.conf import settings
from django import template
import vcf

from main.model_utils import ensure_exists_0775_dir
from main.model_utils import get_dataset_with_type
from main.models import Dataset
from utils import ensure_line_lengths


# Compile this SNPEFF parsing regex only once
# Here we compile the SnpEff-parsing regex that parses a string that looks like
# e.g.:
#   NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aTg/aCg|M239T|386|ygiC||CODING
#        |b3038|1|1)
# Since it's complicated, we build it up in steps:

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
        r'(?P<{:s}>.*)\((?P<{:s}>[^\|]*)',
        r'\|(?P<{:s}>[^\|]*)' * (len(SNPEFF_FIELDS.keys())-4),
        r'\|?(?P<{:s}>[^\|]*)\|?(?P<{:s}>[^\|]*)\)']
        ).format(*SNPEFF_FIELDS.keys()))

MAP_VCF_SOURCE_TOOL_TO_ORIGINAL_VCF_DATASET_TYPE = {
    # TODO: Use constants once circular imports issue is resolved.
    'freebayes': Dataset.TYPE.VCF_FREEBAYES,
    'lumpy': Dataset.TYPE.VCF_LUMPY,
    'pindel': Dataset.TYPE.VCF_PINDEL
}


def build_snpeff(ref_genome):
    """Setup the SnpEff database for ref_genome.

    This function does the following:
        * Sets up the directory structure for SnpEff-related files.
        * Writes a possibly modified Genbank to the location that SnpEff
              expects to find it. A few cleanups are necessary to avoid SnpEff
              quirks.
        * Creates the SnpEff config file for building the database/index.
        * Builds the SnpEff database/index.

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
        print "Skipping SnpEff indexing: No genbank for reference genome %s" % (
                ref_genome.uid)
        return

    # Get the path to the reference genbank, making sure it exists.
    ref_genome_path = get_dataset_with_type(ref_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()
    assert ref_genome_path is not None, "Reference Genbank missing."

    # Create the snpeff directory structure.
    ref_genome.ensure_snpeff_dir()

    # Build a template data dictionary which will be passed to the django
    # template renderer in order to generate the config file.
    templ_data = {}
    templ_data['snpeff_dir'] = ref_genome.get_snpeff_dir()
    templ_data['uid'] = ref_genome.uid
    templ_data['label'] = ref_genome.label

    # The following block does 2 things:
    #    1. Identifies all chromosomes in the Genbank.
    #    2. Ensures that the contained SeqRecord name and ids match, which is
    #       required by SnpEff.
    templ_data['chromosomes'] = []
    new_genbank_seq_records = []
    with open(ref_genome_path) as genbank_fh:
        for seq_record in SeqIO.parse(genbank_fh, 'genbank'):
            # Set the ACCESSION/LOCUS/VERSION to all be the same for this
            # new modified genbank
            seq_record.name = seq_record.id
            new_genbank_seq_records.append(seq_record)

            # Add this record as a chromosome to this ref genome
            # TODO: Do we want to check seqrecords for sane/sanitized names?
            templ_data['chromosomes'].append(seq_record.name)

    templ_data['chromosomes'].append(seq_record.name)
    templ_data['chrs_string'] = ','.join(templ_data['chromosomes'])

    # Write the updated Genbank.
    snpeff_genbank_path = ref_genome.get_snpeff_genbank_file_path()
    SeqIO.write(new_genbank_seq_records, snpeff_genbank_path, 'genbank')

    # Stop-gap fix to ensure line lengths in Genbank to appease SnpEff.
    ensure_line_lengths(ref_genome.get_snpeff_genbank_file_path())

    # Render SnpEff config template.
    render_snpeff_config(templ_data, ref_genome.get_snpeff_config_path())

    # Build snpEff database
    build_snpeff_db(ref_genome.get_snpeff_config_path(), ref_genome.uid)


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

    subprocess.call(snpeff_args)


def get_snpeff_vcf_output_path(alignment_group, vcf_source_tool):
    """Returns the path to SnpEff dir for the given AlignmentGroup and tool.

    Ensures that the intermediate dirs are created so the file can be written.

    Path is of the form:
    /projects/<project_uid>/alignment_groups/vcf/<tool>/snpeff/<tool>_snpeff.vcf
    """
    assert vcf_source_tool in MAP_VCF_SOURCE_TOOL_TO_ORIGINAL_VCF_DATASET_TYPE

    # Root vcf dir.
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)

    # Tool-specific dir.
    tool_dir = os.path.join(vcf_dir, vcf_source_tool)
    ensure_exists_0775_dir(tool_dir)

    # Dir where snpeff data will be written.
    snpeff_vcf_dir = os.path.join(tool_dir, 'snpeff')
    ensure_exists_0775_dir(snpeff_vcf_dir)

    # Snpeff output vcf path.
    vcf_output_filename = os.path.join(
            snpeff_vcf_dir, vcf_source_tool + '_snpeff.vcf')
    return vcf_output_filename


def run_snpeff(alignment_group, vcf_source_tool):
    """Run snpeff on an alignment group after creating a vcf with a snpcaller.

    We only use the alignment type to store the snpeff file.

    Returns the snpeff vcf output filename.
    """
    assert vcf_source_tool in MAP_VCF_SOURCE_TOOL_TO_ORIGINAL_VCF_DATASET_TYPE

    # Get the reference genome uid to get the config path and snpeff genome name
    ref_genome = alignment_group.reference_genome
    ref_genome_uid = alignment_group.reference_genome.uid

    source_vcf_dataset_type = (
            MAP_VCF_SOURCE_TOOL_TO_ORIGINAL_VCF_DATASET_TYPE[vcf_source_tool])
    source_vcf_dataset = get_dataset_with_type(alignment_group,
            type=source_vcf_dataset_type)
    assert source_vcf_dataset is not None
    vcf_input_filename = source_vcf_dataset.get_absolute_location()
    assert os.path.exists(vcf_input_filename)

    # Make sure vcf has at least one record. If not, return.
    with open(vcf_input_filename) as unannotated_fh:
        vcf_reader = vcf.Reader(unannotated_fh)
        try:
            vcf_reader.next()
        except StopIteration:
            # No variants called. No need to do SnpEff.
            return

    # Prepare a directory to put the output file.
    vcf_output_filename = get_snpeff_vcf_output_path(alignment_group,
            vcf_source_tool)

    snpeff_args = [
        'java',
        '-jar', settings.SNPEFF_JAR_PATH,
        'eff',
        '-v',
        '-i',
        'vcf',
        '-o',
        'vcf',
        '-c', ref_genome.get_snpeff_config_path(),
        '-ud', str(settings.SNPEFF_UD_INTERVAL_LENGTH),
        '-formatEff',
        '-q',
        '-noLog',
#        '-t', str(settings.SNPEFF_THREADS),
        ref_genome_uid,
        vcf_input_filename
    ]

    print ' '.join(snpeff_args)

    with open(vcf_output_filename, 'w') as fh_out:
        snpeff_proc = subprocess.Popen(
            snpeff_args,
            stdout=subprocess.PIPE)
        convert_snpeff_info_fields(snpeff_proc.stdout, fh_out)

    return vcf_output_filename


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
    """This function takes a single VCF record and separates the single EFF key
    from snpEFF into multiple separate key-value pairs.

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
    assert hasattr(vcf_record, 'INFO'), 'No INFO attr, not a vcf record.'

    # Dicionary mapping from EFF key to list of values. This will be parsed
    # from the single EFF field output by SnpEff.
    eff_field_lists = defaultdict(list)

    # Check that VCF record has an EFF INFO field
    if 'EFF' in vcf_record.INFO:
        eff_concat_string = vcf_record.INFO['EFF']
    else:
        print >> sys.stderr, ('VCF record at {chrom} {pos} has no '
                'INFO EFF field. Cannot annotate.').format(
                chrom=vcf_record.CHROM, pos=vcf_record.POS)

        vcf_record.INFO['EFF'] = [(
            'ERROR(|||||||||||'
            'SNPEFF_ERROR:NO_EFF_INFO_FIELD|'
            'SNPEFF_ERROR:NO_EFF_INFO_FIELD)')]
        eff_concat_string = vcf_record.INFO['EFF']

    # One EFF group per ALT.
    eff_group_list = []
    for eff in eff_concat_string:
        regex_match = SNPEFF_ALT_RE.match(eff)
        assert regex_match is not None, (
                "Could not parse SnpEff EFF value %s" % eff)
        eff_group_list.append(regex_match)

    eff_fields = list(chain.from_iterable((
            eff.groupdict().items() for eff in eff_group_list)))

    for k, v in eff_fields:
        # mark empty fields as 'bad'
        if v == '':
            v = '.'
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
