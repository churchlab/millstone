"""
Methods for filtering over Variants.
"""

from collections import OrderedDict
import re

from django.db.models import Q
from sympy.core.function import sympify


###############################################################################
# Hard-coded query keys.
# These are hard-coded to the model definitions in main/models.py.
###############################################################################

VARIANT_KEY_MAP = {
    'chromosome': {'type': 'String', 'num': 1},
    'position': {'type': 'Integer', 'num': 1},
}

TYPE_TO_SUPPORTED_OPERATIONS_HARD_CODED = {
        'Float': ['=', '!=', '>=', '<=', '>', '<'],
        'Integer': ['=', '==', '!=', '>=', '<=', '>', '<'],
        'String': ['=', '==', '!='],
        'Boolean': ['=', '==', '!=']
}


###############################################################################
# Query string parsing.
###############################################################################

# Mapping from filter string separator to the django query object
# postfix representation.
# NOTE: Order matters since, for example, we want to try to match '<='
# before we try to match '=' as separator.
DELIM_TO_Q_POSTFIX = OrderedDict()
DELIM_TO_Q_POSTFIX['=='] = ''
DELIM_TO_Q_POSTFIX['<='] = '__lte'
DELIM_TO_Q_POSTFIX['>='] = '__gte'
DELIM_TO_Q_POSTFIX['<'] = '__lt'
DELIM_TO_Q_POSTFIX['>'] = '__gt'
DELIM_TO_Q_POSTFIX['='] = ''


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


def symbol_generator():
    """Generates symbols for manipulation.
    """
    final_symbol_ord = ord('z')
    current_char = 'A'
    while True:
        yield current_char
        next_ord = ord(current_char) + 1
        if next_ord > final_symbol_ord:
            raise StopIteration()
        current_char = chr(next_ord)


EXPRESSION_REGEX = re.compile('(\w+[=><]{1}[=]{0,1}\w+)')

class VariantFilterEvaluator(object):
    """Object that encapsulates methods for evaluating a filter.
    """

    def __init__(self, raw_filter_string, ref_genome):
        self.raw_filter_string = raw_filter_string
        self.clean_filter_string = raw_filter_string.replace(' ', '')
        self.ref_genome = ref_genome
        self._create_symbolic_representation()


    def _create_symbolic_representation(self):
        """Creates a symbolic representation of the query to enable, among
        other things, manipulation with sympy so we can get to disjunctive
        normal form (DNF).
        """
        # Find all the expressions and replace them with symbols.
        self.symbol_to_expression_map = {}

        tokens = EXPRESSION_REGEX.split(self.clean_filter_string)
        symbol_maker = symbol_generator()
        # Iterate through tokens and replace them with successive symbols.
        symbolified_tokens = []
        for token in tokens:
            if EXPRESSION_REGEX.match(token):
                symbol = symbol_maker.next()
                self.symbol_to_expression_map[symbol] = token
                symbolified_tokens.append(symbol)
            else:
                symbolified_tokens.append(token)
        symbolified_string = ''.join(symbolified_tokens)
        self.sympy_representation = sympify(symbolified_string)


    def get_condition_string_for_symbol(self, symbol):
        if isinstance(symbol, str):
            symbol_str = symbol
        else:
            symbol_str = str(symbol)
        return self.symbol_to_expression_map[symbol_str]


class QueryWrapper(object):
    """Object that wraps a simple query, along with metadata depending on the
    target model for the query.
    """

    def __init__(self, raw_query_string):
        self.raw_query_string = raw_query_string

###############################################################################
# Helper methods
###############################################################################

def _get_delim_key_value_triple(raw_string):
    """Attempt to parse a (delim, key, value) triple out of raw_string."""
    # Try the possible delimiters in order until we find one, or fail.
    for raw_delim in DELIM_TO_Q_POSTFIX.iterkeys():
        split_result = raw_string.split(raw_delim)
        delimeter = _clean_delim(raw_delim)
        if len(split_result) == 2:
            key, value = split_result
            for data_map in [VARIANT_KEY_MAP]:
                # Make sure this is a valid key and valid delimeter.
                if key in data_map:
                    specs = data_map[key]
                    if specs['num'] == 1:
                        return tuple([delimeter] + split_result)
                    else:
                        raise ParseError(raw_string, 'Key type not yet supported.')
            # If we got here, the key was not found in any data_map.
            raise ParseError(raw_string, 'Unrecognized filter key.')

    # If we got here, we didn't find a match.
    raise ParseError(raw_string, 'No valid filter delimeter.')


def _clean_delim(raw_delim):
    """Cleans a run delimiter.
    """
    if raw_delim == '=':
        return '=='
    else:
        return raw_delim


def _get_django_q_object_for_triple(delim_key_value_triple):
    """Returns a Django Q object for querying against the SNP model for
    the given key_string.

    Args:
        delim_key_value_triple: A tuple representing a single condition.

    Returns a django Q object.
    """
    assert len(delim_key_value_triple) == 3
    (delim, key, value) = delim_key_value_triple
    postfix = DELIM_TO_Q_POSTFIX[delim]
    return eval('Q(' + key + postfix + '=' + '"' + value + '"' + ')')


###############################################################################
# Main entry method.
###############################################################################

def get_variants_that_pass_filter(filter_string, ref_genome):
    """Takes a complete filter string and returns the variants that pass the
    filter.

    We have a hybrid implementation for field values, where some of the field
    values are actually columns in the respective database tables, while others
    are serialized (pickled) as part of a single key-value string. We do
    filtering over the key-value parts in-memory in the Django server, but on
    a set of objects returned from the SQL database after filtering down by
    the supported DB keys.

    Args:
        filter_string: Query string from the user.
        ref_genome: The ReferenceGenome this is limited to. This is a hard-coded
            parameter of sorts, since at the least we always limit comparisons
            among Variant objects to those that share a ReferenceGenome.

    Returns:
        List of Variant model objects.
    """
    evaluator = VariantFilterEvaluator(filter_string, ref_genome)
    q_obj_list = []
    if not evaluator.sympy_representation.is_Boolean:
        condition_string = evaluator.get_condition_string_for_symbol(
                evaluator.sympy_representation)
        delim_key_value_triple = _get_delim_key_value_triple(condition_string)
        q_obj_list.append(_get_django_q_object_for_triple(
                delim_key_value_triple))
    else:
        for symbol in evaluator.sympy_representation.args:
            condition_string = evaluator.get_condition_string_for_symbol(symbol)
            delim_key_value_triple = _get_delim_key_value_triple(condition_string)
            q_obj_list.append(_get_django_q_object_for_triple(
                    delim_key_value_triple))

    # Combine the Q objects.
    full_q = Q()
    for q_obj in q_obj_list:
        full_q &= q_obj

    # Send the query to the database.
    query_result = ref_genome.variant_set.filter(full_q)
    return query_result
