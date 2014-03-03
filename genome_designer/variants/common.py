"""
Helpers and utils for working with Variants.
"""

from collections import OrderedDict
import re


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
