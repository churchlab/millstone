from collections import OrderedDict
import re

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


SNPEFF_ALT_RE = re.compile(r''.join([
        r'(?P<{:s}>.+)\((?P<{:s}>[^\|]*)',
        r'\|(?P<{:s}>[^\|]*)' * (len(SNPEFF_FIELDS.keys())-4),
        r'\|?(?P<{:s}>[^\|]*)\|?(?P<{:s}>[^\|]*)\)']
        ).format(*SNPEFF_FIELDS.keys()))

s1 = 'missense_variant(MODERATE|MISSENSE|aTg/aCg|p.Met239Thr/c.716T>C|386|ygiC|protein_coding|CODING|b3038|1|C)'
print SNPEFF_ALT_RE.match(s1).groupdict()

s2 = 'stop_lost+splice_region_variant(HIGH|MISSENSE|Taa/Caa|p.Ter387Glnext*?/c.1159T>C|386|ygiC|protein_coding|CODING|b3038|1|C)'
print SNPEFF_ALT_RE.match(s2).groupdict()
