"""
Helpers and utils for working with Variants.
"""

from collections import defaultdict
from collections import OrderedDict
import re

from django.db import connection
from django.db.models import Q

from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantEvidence
from main.models import VariantSet
from main.models import VariantToVariantSet
from variants.filter_eval_result import metadata_default_dict_factory_fn
from variants.filter_eval_result import FilterEvalResult


###############################################################################
# DEBUG (uncomment)
#
# NOTE: Do not commit logging statements unless you know that you want them
#       to be around for a good reason.
###############################################################################

# import logging
# LOGGER = logging.getLogger('debug_logger')
# LOGGER.debug('This adds a logging statement.')


###############################################################################
# Constants
###############################################################################

# Q-object types. This is used to determine whether we need to perform a
# query per sample.
Q_OBJECT_TYPE__GLOBAL = 'global'
Q_OBJECT_TYPE__PER_SAMPLE = 'per_sample'

# SQL key constants
VARIANT_TABLE_KEY__ID = 'id'
VARIANT_TABLE_KEY__SAMPLE = (
        'variantalternate__variantevidence__experiment_sample')

# SQL fields describing a Variant
VARIANT_MODEL_FIELDS = ('id', 'uid', 'type', 'reference_genome_id',
        'chromosome', 'position', 'ref_value')

###############################################################################
# Hard-coded query keys.
# These are hard-coded to the model definitions in main/models.py.
#     * sql_key is relative to Variant model.
###############################################################################

# Keys corresponding to columns in the Variant model.
MELTED_VARIANT_SQL_KEY_MAP = {
    'chromosome': {
        'type': 'String',
        'num': 1,
    },
    'position': {
        'type': 'Integer',
        'num': 1,
    },
    'variant_set_uid': {
        'type': 'String',
        'num': 1,
    },
}

ALL_SQL_KEY_MAP_LIST = [
    MELTED_VARIANT_SQL_KEY_MAP
]

TYPE_TO_SUPPORTED_OPERATIONS = {
        'Float': ['=', '!=', '>=', '<=', '>', '<'],
        'Integer': ['=', '==', '!=', '>=', '<=', '>', '<'],
        'String': ['=', '==', '!='],
        'Boolean': ['=', '==', '!=']
}



################################################################################
# Parsing Regular Expressions
################################################################################

# Mapping from filter string separator to the django query object
# postfix representation.
# NOTE: Order matters since, for example, we want to try to match '<='
# before we try to match '=' as separator.
DELIM_TO_Q_POSTFIX = OrderedDict()
DELIM_TO_Q_POSTFIX['=='] = ''
DELIM_TO_Q_POSTFIX['<='] = '__lte'
DELIM_TO_Q_POSTFIX['>='] = '__gte'
DELIM_TO_Q_POSTFIX['!='] = 'see special handling'
DELIM_TO_Q_POSTFIX['<'] = '__lt'
DELIM_TO_Q_POSTFIX['>'] = '__gt'
DELIM_TO_Q_POSTFIX['='] = ''

# Picks out expressions like, e.g. '(position > 100) in ALL(sample1)'.
CONDITION_PART = '\(.*\)'
SCOPE_TYPE_PART = '(?:ALL|all|ANY|any|ONLY|only)'
SAMPLE_LIST_PART = '\(.*\)'
SAMPLE_SCOPE_REGEX = re.compile(
        '(' + CONDITION_PART + '\s+in\s+' + SCOPE_TYPE_PART +
        SAMPLE_LIST_PART + ')')

# And a named version for grabbing the salient parts.
CONDITION_PART_NAMED = '\((?P<condition>.*)\)'
SCOPE_TYPE_PART_NAMED = '(?P<scope_type>ALL|all|ANY|any|ONLY|only)'
SAMPLE_LIST_PART_NAMED = '\((?P<samples>.*)\)'
SAMPLE_SCOPE_REGEX_NAMED = re.compile(
        '(' + CONDITION_PART_NAMED + '\s+in\s+' + SCOPE_TYPE_PART_NAMED +
        SAMPLE_LIST_PART_NAMED + ')')

# Recognizes a pattern of the form 'key op value'.
EXPRESSION_REGEX = re.compile('(\w+\s*[=><!]{1}[=]{0,1}\s*\w+)')

# Recognizes statements about set.
SET_REGEX = re.compile('((?:NOT_){0,1}IN_SET\([\w]+\))')
SET_REGEX_NAMED = re.compile('((?P<maybe_not>NOT_){0,1}IN_SET\((?P<sets>[\w]+)\))')

GENE_REGEX = re.compile('(IN_GENE\([\w]+\))')
GENE_REGEX_NAMED = re.compile('(IN_GENE\((?P<gene>[\w]+)\))')


###############################################################################
# Helper objects.
###############################################################################

class ParseError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg

    def __str__(self):
        return ('Parse Error while parsing: ' + str(self.expr) + '. ' +
                str(self.msg))

def SymbolGenerator():
    """Generator that yields characters in alphabetical order, starting with
    'A' and continuing on to 'Z', followed by 'a' and ending with 'z'.

    NOTE: The current implementation runs out of symbols after 26 * 2 symbols.
    """
    final_symbol_ord = ord('z')
    current_char = 'A'
    while True:
        yield current_char
        next_ord = ord(current_char) + 1
        if next_ord > final_symbol_ord:
            raise StopIteration()
        current_char = chr(next_ord)


class SqlReadySymbol(object):
    """An object that packages a representation of a symbol that is supported
    for SQL queries into an object that provides the Django Q object as well
    as the raw (delim, key, value) representation for in-memory processing.
    """

    def __init__(self, semantic_type, q_obj, delim_key_value_triple=None):
        self.semantic_type = semantic_type
        self.q_obj = q_obj

        # TODO: This is not consistently useful for all semantic types.
        # Figure out more robust common solution.
        self.delim_key_value_triple = delim_key_value_triple


################################################################################
# Utility Methods
################################################################################

def get_delim_key_value_triple(raw_string, all_key_map):
    """Attempt to parse a (delim, key, value) triple out of raw_string."""
    # Remove spaces from the string.
    raw_string = raw_string.replace(' ', '')

    # Try the possible delimiters in order until we find one, or fail.
    for raw_delim in DELIM_TO_Q_POSTFIX.iterkeys():
        split_result = raw_string.split(raw_delim)
        delimeter = _clean_delim(raw_delim)
        if len(split_result) == 2:
            key, value = split_result
            for data_map in all_key_map:
                # Make sure this is a valid key and valid delimeter.
                if key in data_map:
                    specs = data_map[key]
                    if specs['num'] in (-1,1):
                        return tuple([delimeter] + split_result)
                    else:
                        raise ParseError(raw_string, 
                                'Key type {:d} not yet supported.'.format(
                                specs['num']))
            # If we got here, the key was not found in any data_map.
            raise ParseError(raw_string, 'Unrecognized filter key: %s' % key)

    # If we got here, we didn't find a match.
    raise ParseError(raw_string, 'No valid filter delimeter.')


def convert_delim_key_value_triple_to_expr(triple):
    (delim, key, value) = triple
    # Make '==' SQL-friendly.
    if delim == '==':
        delim = '='
    return (key + delim + '%s', value)


def _clean_delim(raw_delim):
    """Cleans a run delimiter.
    """
    if raw_delim == '=':
        return '=='
    else:
        return raw_delim


def get_all_key_map(ref_genome):
    return ALL_SQL_KEY_MAP_LIST + ref_genome.variant_key_map.values()


def extract_filter_keys(filter_expr, ref_genome):
    """Returns a list of keys that we are filtering over in the expression.
    """
    filter_keys = []
    for expr in EXPRESSION_REGEX.finditer(filter_expr):
        (delim, key, value) = get_delim_key_value_triple(expr.group(),
                get_all_key_map(ref_genome))
        filter_keys.append(key)
    return filter_keys


def get_sample_id_set_for_variant(variant):
    """Returns the set of all ExperimentSamples ids for which there exists a
    relationship to the Variant
    """
    return set([ve.experiment_sample_id for ve in
            VariantEvidence.objects.filter(
                    variant_caller_common_data__variant=variant)])


def get_django_q_object_for_triple(delim_key_value_triple):
    """Returns a Django Q object for querying against the SNP model for
    the given key_string.

    Args:
        delim_key_value_triple: A tuple representing a single condition.

    Returns a django Q object.
    """
    assert len(delim_key_value_triple) == 3
    (delim, key, value) = delim_key_value_triple

    # Default type for this Q object.
    q_object_type = Q_OBJECT_TYPE__GLOBAL

    # Figure out proper model to add a model prefix in the Q object, e.g.:
    #     Q(variantalternate__alt_value=A)
    model_prefix = ''
    if key in VARIANT_ALTERNATE_SQL_KEY_MAP:
        model_prefix = 'variantalternate__'
        q_object_type = Q_OBJECT_TYPE__PER_SAMPLE

    # Special handling for != delim.
    if delim == '!=':
        postfix = ''
        maybe_not_prefix = '~'
    else:
        postfix = DELIM_TO_Q_POSTFIX[delim]
        maybe_not_prefix = ''

    eval_string = (maybe_not_prefix + 'Q(' + model_prefix + key + postfix +
            '=' + '"' + value + '"' + ')')
    q_obj = eval(eval_string)
    return SqlReadySymbol(q_object_type, q_obj,
            delim_key_value_triple=delim_key_value_triple)


def get_django_q_object_for_gene_restrict(gene_restrict_string, ref_genome):
    """Returns a tuple (Q_OBJECT_TYPE, Q object to limit results to the gene).
    """
    match = GENE_REGEX_NAMED.match(gene_restrict_string)
    gene_label = match.group('gene')

    # Restrict to variants whose position fall between start and end of
    # gene.

    # First attempt, look up the Gene and get positions.
    gene_region = Region.objects.get(
            type=Region.TYPE.GENE,
            reference_genome=ref_genome,
            label=gene_label)

    # Assume gene only has one interval.
    gene_interval = gene_region.regioninterval_set.all()[0]

    # Return a Q object bounding Variants by this position.
    q_obj = (Q(position__gte=gene_interval.start) &
            Q(position__lt=gene_interval.end))
    return SqlReadySymbol(Q_OBJECT_TYPE__GLOBAL, q_obj)


def evaluate_condition_in_triple(data_map, type_map, triple, idx=None):
    """Evaluates a condition.

    idx arg loops through all possible items by list index if the data type is a
    list of values (i.e. in the per alternate case - if the data_type map 'spec'
    field is '-1', corresponding to a Number='A' in the vcf). If it is empty,
    then evaluate the data type as one value.

    Idx field calls are  recursive calls from within the function to iterate
    through the list. If any values are true, then the condition returns true.
    """

    (delim, key, value) = triple

    # If this is an INFO field (common_data) and it is a per-alternate field
    # (Number = 'A' in vcf, 'num' == -1 in pyvcf), then match if any of the
    # values is correct. This recursively calls the function with the extra idx
    # field.
    if idx is None and 'num' in type_map[key] and type_map[key]['num'] == -1:
        evaluations = []
        for recurse_idx in range(len(data_map[key])):
            evaluations.append(evaluate_condition_in_triple(
                    data_map,
                    type_map,
                    triple, 
                    idx=recurse_idx))
        return any(evaluations)
    else:
        cast_type_string = type_map[key]['type']
        if cast_type_string == 'Boolean':
            return _evaluate_boolean_condition(data_map, key, value, idx)
        else:
            casted_value = _cast_value_to_type(value, cast_type_string)
            if idx is not None:
                evaled_list = data_map[key]
                return eval('evaled_list[idx] ' + delim + ' casted_value')
            else:
                return eval('data_map[key] ' + delim + ' casted_value')

def _evaluate_boolean_condition(data_dict, key, value, idx=None):
    """Evaluates a boolean condition.
    """
    VALID_BOOLEAN_TRUE_VALUES = ['True', 'true', 'T', 't']
    VALID_BOOLEAN_FALSE_VALUES = ['False', 'false', 'F', 'f']

    #if data_dict[key] is a dictionary
    if idx is not None:
        init_result = eval(data_dict[key])[idx]
    else:
        init_result = data_dict[key]
    if value in VALID_BOOLEAN_TRUE_VALUES:
        return init_result
    elif value in VALID_BOOLEAN_FALSE_VALUES:
        return not init_result
    else:
        raise ParseError(value, 'Invalid boolean value, use True or False')


def _cast_value_to_type(value, cast_type_string):
    """Return the value casted to the type specified by the cast_type_string,
    as defined in the type maps in the ReferenceGenome object's variant_key_map field
    """
    if cast_type_string == 'Integer':
        return int(value)
    elif cast_type_string == 'Float':
        return float(value)
    elif cast_type_string == 'String':
        return str(value)
    else:
        raise Exception("Unsupported type " + cast_type_string)


def assert_delim_for_key(type_map, delim, key):
    """Asserts that the delimiter can be evaluated for the type comparison
    specified by the key. Raises a ParseError if not.
    """
    data_type = type_map[key]['type']
    if not delim in TYPE_TO_SUPPORTED_OPERATIONS[data_type]:
        raise ParseError(str(key) + str(delim),
                'Invalid delim for type indicated by key.')


def get_variant_table_column_for_sql_key(sql_key):
    """Returns the name of the column in the Variant table for the sql key.

    Raises:
        AssertionError if key is not recognized.
    """
    for sql_key_map in ALL_SQL_KEY_MAP_LIST:
        if sql_key in sql_key_map:
            return sql_key_map[sql_key]['variant_table_col']
    raise AssertionError("Invalid sql_key: %s" % sql_key)


def dictfetchall(cursor):
    """Returns all rows from a cursor as a dict.

    Source: https://docs.djangoproject.com/en/dev/topics/db/sql/
    """
    desc = cursor.description
    return [
        dict(zip([col[0] for col in desc], row))
        for row in cursor.fetchall()
    ]


class HashableVariantDict(object):
    """Hashable version of Variant Dict.

    The Variant id and Sample id uniquely identify the row.
    """

    def __init__(self, obj_dict):
        self.obj_dict = obj_dict

    def __hash__(self):
        return hash((self.obj_dict['id'],
            self.obj_dict['experiment_sample_id']))

    def __getitem__(self, key):
        return self.obj_dict[key]

    def __getattr__(self, attr):
        """Override.

        Since we changed the code from using python objects and instead
        dictionaries, there are plenty of places using 'dot' attribute access
        so rather than looking for all of those places I'm going to try making
        this small hack and see how it works out.
        """
        try:
            return self.obj_dict[attr]
        except KeyError:
            raise AttributeError


def hashablefetchall(cursor):
    """Returns all rows from a cursor as hashable object.
    """
    desc = cursor.description
    return [
        HashableVariantDict(dict(zip([col[0] for col in desc], row)))
        for row in cursor.fetchall()
    ]


def eval_variant_set_filter_expr(set_restrict_string, ref_genome):
    """Returns a FilterEvalResult containing Variants and per-sample metadata
    that satisfy the passed in set filter.

    Args:
        set_restrict_string: Valid expression for a set.

    Returns:
        A FilterEvalResult object.
    """
    return _eval_variant_set_filter_expr__brute_force(
            set_restrict_string, ref_genome)


def _eval_variant_set_filter_expr__brute_force(
        set_restrict_string, ref_genome):
    """Implementation that uses brute force on negative case.

    This will help us come up with a bunch of tests for when we implement the
    fully optimized version.

    Args:
        set_restrict_string: The filter token (e.g. 'IN_SET(abcd1234)').
        ref_genome: ReferenceGenome this search is relative to.

    Returns:
        A FilterEvalResult object.

    Raises:
        DoesNotExist: If VariantSet with request uid doesn't exist.
    """
    # Extract the relevant part of the query.
    match = SET_REGEX_NAMED.match(set_restrict_string)
    variant_set_uid = match.group('sets')
    assert len(variant_set_uid) > 0, (
            "No actual set provided in set filter.")

    # Look up the VariantSet.
    variant_set = VariantSet.objects.get(
            uid=variant_set_uid,
            reference_genome=ref_genome)

    # Determine whether this is a NOT query.
    is_negative_query = bool(match.group('maybe_not'))

    # If forward-case, use optimized implementation.
    if is_negative_query:
        return _eval_variant_set_filter_expr__brute_force__negative(
                variant_set.id, ref_genome)
    else:
        return _eval_variant_set_filter_expr__brute_force__positive(
                variant_set.id, ref_genome)


def _eval_variant_set_filter_expr__brute_force__positive(variant_set_id,
        ref_genome):
    """Positive look up by brute force.

    NOTE: We eventually want to hard-code the SQL for this as started in
    _eval_variant_set_filter_expr__optimized__positive() but this broke
    when we switched to Postgres since it has a slightly different syntax
    from Sqlite.
    """
    # Store the results in these data structures.
    passing_variant_set = set()
    variant_id_to_metadata_dict = defaultdict(metadata_default_dict_factory_fn)

    for variant in Variant.objects.filter(reference_genome=ref_genome):
        matching_vtvs = VariantToVariantSet.objects.filter(
                variant=variant,
                variant_set_id=variant_set_id)
        if len(matching_vtvs) > 0:
            assert len(matching_vtvs) == 1, (
                    "Multiple VariantToVariantSet objects found for "
                    "Variant.uid=%s, VariantSet.uid=%s" % (
                            variant.uid, variant_set.uid))
            vtvs = matching_vtvs[0]
            passing_variant_set.add(variant)
            passing_sample_ids = set([sample.id for sample in
                    vtvs.sample_variant_set_association.all()])
            variant_id_to_metadata_dict[variant.id][
                    'passing_sample_ids'] = passing_sample_ids

    return FilterEvalResult(passing_variant_set, variant_id_to_metadata_dict)


def _eval_variant_set_filter_expr__brute_force__negative(variant_set_id,
        ref_genome):
    """Negative look up by brute force.

    NOTE: We eventually want to hard-code the SQL for this as started in
    _eval_variant_set_filter_expr__optimized_full_not_working().
    """
    # Store the results in these data structures.
    passing_variant_set = set()
    variant_id_to_metadata_dict = defaultdict(metadata_default_dict_factory_fn)

    # Otherwise, do brute force solution that
    for variant in Variant.objects.filter(reference_genome=ref_genome):
        # Add an entry only if at least one passes the filter (of not
        # being in this set).
        passing_sample_ids = set()
        has_at_least_one_evidence_obj = False
        for common_data_obj in variant.variantcallercommondata_set.all():
            for variant_evidence in common_data_obj.variantevidence_set.all():
                has_at_least_one_evidence_obj = True
                experiment_sample = variant_evidence.experiment_sample
                matching_vtvs = VariantToVariantSet.objects.filter(
                        variant=variant,
                        variant_set_id=variant_set_id,
                        sample_variant_set_association__id=experiment_sample.id)
                if not len(matching_vtvs):
                    passing_sample_ids.add(experiment_sample.id)
        if passing_sample_ids:
            passing_variant_set.add(variant)
            variant_id_to_metadata_dict[variant.id][
                    'passing_sample_ids'] = passing_sample_ids
        elif not has_at_least_one_evidence_obj:
            # Another way for this Variant to pass is if no VariantEvidence
            # objects exist for the Variant in which case we do a simpler
            # query.
            matching_vtvs = VariantToVariantSet.objects.filter(
                    variant=variant, variant_set_id=variant_set_id)
            if not len(matching_vtvs):
                passing_variant_set.add(variant)

    return FilterEvalResult(passing_variant_set, variant_id_to_metadata_dict)


def _eval_variant_set_filter_expr__optimized__positive(variant_set_id,
        ref_genome):
    """Supports positive queries in a more optimal fashion.

    NOTE: Temporarily unused as this broke when we switched to Postgresql.
    """
    # We need to perform a custom SQL statement to get all the INNER JOINs
    # right. As far as we can tell, there is not a clear way to do this using
    # Django's ORM.
    cursor = connection.cursor()

    sql_statement = (
        'SELECT '
            # SELECT Variant fields
            '"main_variant"."id", "main_variant"."uid", '
            '"main_variant"."type", "main_variant"."reference_genome_id", '
            '"main_variant"."chromosome", "main_variant"."position", '
            '"main_variant"."ref_value", '

            # SELECT ExperimentSample.id
            '"main_varianttovariantset_sample_variant_set_association"."experimentsample_id" AS "sample_id" '

        # FROM a bunch of JOINed tables
        'FROM "main_variant"'
            'INNER JOIN "main_varianttovariantset" '
                'ON ("main_variant"."id" = "main_varianttovariantset"."variant_id") '
            'INNER JOIN "main_variantset" '
                'ON ("main_varianttovariantset"."variant_set_id" = "main_variantset"."id") '
            'LEFT OUTER JOIN "main_varianttovariantset_sample_variant_set_association" '
                'ON ("main_varianttovariantset"."id" = "main_varianttovariantset_sample_variant_set_association"."varianttovariantset_id") '

        # WHERE VariantSet.uid is what is in the filter expr
        'WHERE ("main_variant"."reference_genome_id" = %s AND '
            '"main_variantset"."id" = "%s")'

        % (ref_genome.id, variant_set_id)
    )

    data_dict_list = dictfetchall(cursor.execute(sql_statement))

    # Gather the Variant objects from the data.
    passing_variant_set = set()
    variant_ids_visited = set()
    for row in data_dict_list:
        if row['id'] in variant_ids_visited:
            continue
        passing_variant_set.add(Variant(**get_subdict(row, VARIANT_MODEL_FIELDS)))
        variant_ids_visited.add(row['id'])

    # Build the metadata dictionary.
    variant_id_to_metadata_dict = defaultdict(metadata_default_dict_factory_fn)
    for row in data_dict_list:
        if row['sample_id']:
            variant_id_to_metadata_dict[row['id']]['passing_sample_ids'].add(
                row['sample_id'])

    return FilterEvalResult(passing_variant_set, variant_id_to_metadata_dict)


def _eval_variant_set_filter_expr__optimized_full_not_working(
        set_restrict_string, ref_genome):
    """Optimized SQL implementation.

    TODO: Finish.
    """
    # Extract the relevant part of the query.
    match = SET_REGEX_NAMED.match(set_restrict_string)
    variant_set_uid = match.group('sets')
    assert len(variant_set_uid) > 0, (
            "No actual set provided in set filter.")

    # Get and validate the VariantSet.
    variant_set = VariantSet.objects.get(
            uid=variant_set_uid,
            reference_genome=ref_genome)

    # We need to perform a custom SQL statement to get all the INNER JOINs
    # right. As far as we can tell, there is not a clear way to do this using
    # Django's ORM.
    cursor = connection.cursor()

    if not match.group('maybe_not'):
        return _eval_variant_set_filter_expr__optimized__positive(
                variant_set.id, ref_genome)
    else:
        # SELECT Variant.id that are in the VariantSet with requested uid.
        inner_select_statement = (
            'SELECT '
                '"main_varianttovariantset"."id" AS "vtvs_id", '
                '"main_varianttovariantset"."variant_id" AS "variant_id", '
                '"main_variantset"."uid" AS "variant_set_uid", '
                '"main_varianttovariantset_sample_variant_set_association"."experimentsample_id" AS "sample_id" '


            'FROM "main_varianttovariantset" '
                'INNER JOIN "main_variantset" '
                    'ON ("main_varianttovariantset"."variant_set_id" = "main_variantset"."id") '
                'LEFT OUTER JOIN "main_varianttovariantset_sample_variant_set_association" '
                    'ON ("main_varianttovariantset"."id" = "main_varianttovariantset_sample_variant_set_association"."varianttovariantset_id") '

            'WHERE "main_variantset"."uid" = "%s"'

            % variant_set_uid
        )

        # # DEBUG
        data_dict_list = dictfetchall(cursor.execute(inner_select_statement))
        print 'INNER_SELECT', len(data_dict_list)
        for row in data_dict_list:
            print row
        # assert False, "STOP"

        sql_statement = (
            'SELECT '
                # SELECT Variant fields
                '"main_variant"."id", "main_variant"."uid", '
                '"main_variant"."type", "main_variant"."reference_genome_id", '
                '"main_variant"."chromosome", "main_variant"."position", '
                '"main_variant"."ref_value", '

                # SELECT ExperimentSample.id
                'U1."sample_id" AS "sample_id" '

            # FROM a bunch of JOINed tables
            'FROM "main_variant" '
                'LEFT OUTER JOIN (' + inner_select_statement + ') U1 '
                    'ON ("main_variant"."id" = U1."variant_id") '
                # 'LEFT OUTER JOIN "main_varianttovariantset_sample_variant_set_association" '
                #     'ON (U1."vtvs_id" = "main_varianttovariantset_sample_variant_set_association"."varianttovariantset_id") '


            'WHERE U1."variant_set_uid" IS NULL'
        )

        # DEBUG
        data_dict_list = dictfetchall(cursor.execute(sql_statement))
        print 'NUM_RESULTS', len(data_dict_list)
        for row in data_dict_list:
            print row

    data_dict_list = dictfetchall(cursor.execute(sql_statement))

    # Gather the Variant objects from the data.
    passing_variant_set = set()
    variant_ids_visited = set()
    for row in data_dict_list:
        if row['id'] in variant_ids_visited:
            continue
        passing_variant_set.add(Variant(**get_subdict(row, VARIANT_MODEL_FIELDS)))
        variant_ids_visited.add(row['id'])

    # Build the metadata dictionary.
    variant_id_to_metadata_dict = defaultdict(metadata_default_dict_factory_fn)
    for row in data_dict_list:
        if row['sample_id']:
            variant_id_to_metadata_dict[row['id']]['passing_sample_ids'].add(
                row['sample_id'])

    return FilterEvalResult(passing_variant_set, variant_id_to_metadata_dict)


def get_subdict(big_dict, fields):
    """Returns a dict that is the subset of the dict with only the specified
    fields.

    Args:
        big_dict: Dictionary with keys that are superset of fields.
        fields: Iterable of field keys.

    Raises:
        TypeError if field missing from big_dict.
    """
    return {k: big_dict[k] for k in fields}


def create_initial_filter_eval_result_object(variant_query_set):
    """Creates a FilterEvalResult object containing the variants
    and all samples associated with those variants.
    """
    variant_id_to_metadata_dict = defaultdict(metadata_default_dict_factory_fn)

    query_result = _get_variant_id_sample_id_tuple_list(variant_query_set)

    for variant_id, sample_id in query_result:
        variant_id_to_metadata_dict[variant_id]['passing_sample_ids'].add(
                sample_id)

    return FilterEvalResult(set(variant_query_set), variant_id_to_metadata_dict)


def _get_variant_id_sample_id_tuple_list(variant_query_set):
    """Returns list of two-tuples (variant_id, sample_id).

    NOTE: Transient implementation that ignores variant_query_set.
    """
    cursor = connection.cursor()
    sql_statement = (
            'SELECT id, experiment_sample_id FROM materialized_melted_variant')
    cursor.execute(sql_statement)
    return cursor.fetchall()
