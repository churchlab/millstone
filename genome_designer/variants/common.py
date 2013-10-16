"""
Helpers and utils for working with Variants.
"""

from collections import OrderedDict
import re

from django.db.models import Q

from main.models import VariantEvidence


###############################################################################
# Constants
###############################################################################

# Q-object types. This is used to determine whether we need to perform a
# query per sample.
Q_OBJECT_TYPE__GLOBAL = 'global'
Q_OBJECT_TYPE__PER_SAMPLE = 'per_sample'


###############################################################################
# Hard-coded query keys.
# These are hard-coded to the model definitions in main/models.py.
###############################################################################

# Keys corresponding to columns in the Variant model.
VARIANT_SQL_KEY_MAP = {
    'chromosome': {'type': 'String', 'num': 1},
    'position': {'type': 'Integer', 'num': 1},
}

# Keys corresponding to columns in the Variant model.
VARIANT_ALTERNATE_SQL_KEY_MAP = {
    'alt_value': {'type': 'String', 'num': 1},
}

# Keys corresponding to columns in the VariantCallerCommonData model.
VARIANT_CALLER_COMMON_DATA_SQL_KEY_MAP = {
}

# Keys corresponding to columns in the VariantEvidence model.
VARIANT_EVIDENCE_SQL_KEY_MAP = {
}

ALL_SQL_KEY_MAP_LIST = [
    VARIANT_SQL_KEY_MAP,
    VARIANT_ALTERNATE_SQL_KEY_MAP,
    VARIANT_CALLER_COMMON_DATA_SQL_KEY_MAP,
    VARIANT_EVIDENCE_SQL_KEY_MAP,
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
    return set([ve.experiment_sample.id for ve in
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
    return (q_object_type, q_obj)


def get_django_q_object_for_set_restrict(set_restrict_string):
    """Returns a tuple (Q_OBJECT_TYPE, Q object limiting reuslts to a set).
    """
    match = SET_REGEX_NAMED.match(set_restrict_string)
    variant_set_uid = match.group('sets')
    assert len(variant_set_uid) > 0, (
            "No actual set provided in set filter.")

    q_obj = Q(varianttovariantset__variant_set__uid=variant_set_uid)

    # Maybe negate the result.
    if match.group('maybe_not'):
        q_obj = ~q_obj

    return (Q_OBJECT_TYPE__GLOBAL, q_obj)


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
    return (Q_OBJECT_TYPE__GLOBAL, q_obj)


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
