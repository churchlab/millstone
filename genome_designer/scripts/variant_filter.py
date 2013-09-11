"""
Methods for filtering over Variants.

The fields that we allow querying against vary across different models,
and for each model, the field may either be a column in the SQL database, a key
in the catch-all serialized key-value dictionary.

A naive implementation would pull all Variant objects from the SQl database
(limited to the ReferenceGenome provided) into memory, and then iterate over
all these objects for each condition in the query. However, we could reduce
the number of objects pulled into memory by extracting the parts of the query
that correspond to columns in the SQL table. The strategy for doing this is
to convert the query into disjunctive normal form, an OR of ANDs, and then
evaluate each conjunction (AND clause) separately. Finally, take the union
of the result of each evaluated conjunction.

We can use sympy, a python symbolic mathematics library, to help convert the
user's query into DNF. The main idea here is to take the raw query, replace
conditions with a symbol, convert to DNF using sympy, and then evaluate each
conjunctive clause as described above.

Other implementation nuances to note:
    * Since we are comparing across multiple tables that are inter-related,
      remember to use Django's select_related() method where appropriate, in
      order to limit the number of independent hits of the SQL db.
"""

from collections import OrderedDict
import re

from django.db.models import Q
from sympy.logic import boolalg


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

# TODO: Generate these from the vcf dataset(s) associated with the reference
# genome we are querying against.
from snp_filter_key_map import VARIANT_CALLER_COMMON_MAP
from snp_filter_key_map import VARIANT_EVIDENCE_MAP

ALL_KEY_MAP_LIST = ALL_SQL_KEY_MAP_LIST + [
    VARIANT_CALLER_COMMON_MAP,
    VARIANT_EVIDENCE_MAP,
]

TYPE_TO_SUPPORTED_OPERATIONS = {
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

        # Generator object that provides symbols in alphabetical order.
        symbol_maker = symbol_generator()

        # Iterate through tokens and replace them with successive symbols.
        tokens = EXPRESSION_REGEX.split(self.clean_filter_string)
        symbolified_tokens = []
        for token in tokens:
            if EXPRESSION_REGEX.match(token):
                symbol = symbol_maker.next()
                self.symbol_to_expression_map[symbol] = token
                symbolified_tokens.append(symbol)
            else:
                symbolified_tokens.append(token)
        symbolified_string = ''.join(symbolified_tokens)
        self.sympy_representation = boolalg.to_dnf(symbolified_string)


    def get_condition_string_for_symbol(self, symbol):
        """Returns the condition string that the symbol replaced.

        Args:
            symbol: A sympy.core.symbol.Symbol or string.

        Returns:
            A string representing a query condition.
        """
        if isinstance(symbol, str):
            symbol_str = symbol
        else:
            symbol_str = str(symbol)
        return self.symbol_to_expression_map[symbol_str]


    def evaluate(self):
        """Evaluates the conditions in the filter string.

        Returns:
            Set of Variants that pass the filter query.
        """
        return self.evaluate_disjunction(self.sympy_representation)


    def evaluate_disjunction(self, disjunction_clause):
        """Evaluate a disjunction (OR) clause.

        Args:
            disjunction_clause: sympy.logic.boolalg.Or clause.

        Returns:
            Set of variants that satisfy the clause.
        """
        # The disjunction may contain a single condition.
        if not isinstance(disjunction_clause, boolalg.Or):
            return set(self.evaluate_conjunction(disjunction_clause))
        else:
            result = set()
            for conjunction_clause in disjunction_clause.args:
                result |= set(self.evaluate_conjunction(conjunction_clause))
            return result


    def evaluate_conjunction(self, conjunction_clause):
        """Evaluate a conjunction condition. That is all symbols are ANDed.
        """
        # Iterate through the conditions corresponding to the symbols in the
        # clause and either create Q objects out of them, relative to the
        # Variant model, or save them as key-value models to evaluate in memory
        # after the SQL fetch. Order doesn't matter since this is a conjunction
        # clause.

        def _single_symbol_helper(symbol, q_part, remaining_triples):
            """Helper method for evaluating a single symbol.
            """
            q_or_triple = self.handle_single_symbol(symbol)
            if isinstance(q_or_triple, Q):
                q_part &= q_or_triple
            else:
                remaining_triples.append(q_or_triple)
            return (q_part, remaining_triples)

        q_part = Q()
        remaining_triples = []

        # The conjunction may contain a single condition.
        if not conjunction_clause.is_Boolean:
            (q_part, remaining_triples) = _single_symbol_helper(
                    conjunction_clause, q_part, remaining_triples)
        else:
            for symbol in conjunction_clause.args:
                (q_part, remaining_triples) = _single_symbol_helper(
                        symbol, q_part, remaining_triples)

        # Now perform the SQL query to pull the models that pass the Q part
        # into memory.
        variant_list = self.ref_genome.variant_set.filter(q_part)

        # NOTE: Getting this correct first. After that, optimize lookups.
        for triple in remaining_triples:
            (delim, key, value) = triple
            passing_variant_list = []
            for variant in variant_list:
                if key in VARIANT_CALLER_COMMON_MAP:
                    _assert_delim_for_key(VARIANT_CALLER_COMMON_MAP, delim, key)
                    all_common_data_obj = (
                            variant.variantcallercommondata_set.all())
                    # TODO: Figure out semantics of having more than one common
                    # data object.
                    for common_data_obj in all_common_data_obj:
                        data_dict = common_data_obj.as_dict()
                        passing = _evaluate_condition_in_triple(
                                data_dict, VARIANT_CALLER_COMMON_MAP, triple)
                        if passing:
                            passing_variant_list.append(variant)
                            break
                else:
                    raise ParseError(key, 'Unrecognized filter key.')

                for common_data_obj in variant.variantcallercommondata_set.all():
                    data_dict = common_data_obj.as_dict()
                    data_dict

            # Since we are inside of a conjunction, we only need to check
            # the remaining variants on the next iteration.
            variant_list = passing_variant_list

        return variant_list


    def handle_single_symbol(self, symbol):
        """Returns one of:
             * A Django Q object if the symbol represents a condition that
                can be evaluated against the SQL database.
             * A triple of delim, key, value if the condition must be evaluated
                in-memory.
        """
        condition_string = self.get_condition_string_for_symbol(symbol)
        (delim, key, value) = _get_delim_key_value_triple(condition_string)
        for key_map in ALL_SQL_KEY_MAP_LIST:
            if key in key_map:
                return _get_django_q_object_for_triple((delim, key, value))
            else:
                return (delim, key, value)


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
            for data_map in ALL_KEY_MAP_LIST:
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


def _evaluate_condition_in_triple(data_map, type_map, triple):
    """Evaluates a condition.
    """
    (delim, key, value) = triple
    cast_type_string = type_map[key]['type']
    if cast_type_string == 'Boolean':
        return _evaluate_boolean_condition(data_map, key, value)
    else:
        casted_value = _cast_value_to_type(value, cast_type_string)
        return eval('data_map[key] ' + delim + ' casted_value')


def _evaluate_boolean_condition(data_dict, key, value):
    """Evaluates a boolean condition.
    """
    VALID_BOOLEAN_TRUE_VALUES = ['True', 'true', 'T', 't']
    VALID_BOOLEAN_FALSE_VALUES = ['False', 'false', 'F', 'f']
    init_result = data_dict[key]
    if value in VALID_BOOLEAN_TRUE_VALUES:
        return init_result
    elif value in VALID_BOOLEAN_FALSE_VALUES:
        return not init_result
    else:
        raise ParseError(value, 'Invalid boolean value, use True or False')


def _cast_value_to_type(value, cast_type_string):
    """Return the value casted to the type specified by the cast_type_string,
    as defined in snp_filter_key_map.py
    """
    if cast_type_string == 'Integer':
        return int(value)
    elif cast_type_string == 'Float':
        return float(value)
    elif cast_type_string == 'String':
        return str(value)
    else:
        raise Exception("Unsupported type " + cast_type_string)


def _assert_delim_for_key(type_map, delim, key):
    """Asserts that the delimiter can be evaluated for the type comparison
    specified by the key. Raises a ParseError if not.
    """
    data_type = type_map[key]['type']
    if not delim in TYPE_TO_SUPPORTED_OPERATIONS[data_type]:
        raise ParseError(str(key) + str(delim),
                'Invalid delim for type indicated by key.')


###############################################################################
# Main client method.
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
    return evaluator.evaluate()
