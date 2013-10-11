"""
Utils for working with Variants.
"""

from collections import OrderedDict

import re


###############################################################################
# Hard-coded query keys.
# These are hard-coded to the model definitions in main/models.py.
###############################################################################

# Keys corresponding to columns in the Variant model.
VARIANT_SQL_KEY_MAP = {
    'chromosome': {'type': 'String', 'num': 1},
    'position': {'type': 'Integer', 'num': 1},
}

# Keys corresponding to columns in the VariantCallerCommonData model.
VARIANT_CALLER_COMMON_DATA_SQL_KEY_MAP = {
}

# Keys corresponding to columns in the VariantEvidence model.
VARIANT_EVIDENCE_SQL_KEY_MAP = {
}

ALL_SQL_KEY_MAP_LIST = [
    VARIANT_SQL_KEY_MAP,
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
