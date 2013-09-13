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

from main.models import ExperimentSample
from main.models import VariantEvidence


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
DELIM_TO_Q_POSTFIX['!='] = 'see special handling'
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


class FilterEvalResult(object):
    """Wraps the result of evaluating a filter condition.

    Provides utility methods for combining results.
    """

    def __init__(self, variant_set, variant_id_to_metadata_dict):
        if not isinstance(variant_set, set):
            variant_set = set(variant_set)
        self.variant_set = variant_set
        self.variant_id_to_metadata_dict = variant_id_to_metadata_dict

    def __or__(self, other):
        return self.combine(other, '|')

    def __and__(self, other):
        return self.combine(other, '&')

    def combine(self, other, op_string):
        """Method that returns a new FilterEvalResult that is the combination
        of this one and other.

        Args:
            other: The FilterEvalResult to combine with.
            op_string: Either '&' or '|'.

        Returns:
            A new FilterEvalResult object.
        """
        assert isinstance(other, FilterEvalResult)
        assert op_string in ['&', '|']

        # Combine the Variant sets.
        if op_string == '&':
            new_variant_set = self.variant_set & other.variant_set
        elif op_string == '|':
            new_variant_set = self.variant_set | other.variant_set
        else:
            raise AssertionError("Unsupported op: %s" % op_string)

        # Build up the new metadata map.
        new_variant_id_to_metadata_dict = {}
        for variant in new_variant_set:
            merged_filter_metadata = {}

            self_filter_metadata = self.variant_id_to_metadata_dict.get(
                    variant.id, {})
            other_filter_metadata = other.variant_id_to_metadata_dict.get(
                    variant.id, {})

            # Merge passing sample ids.
            self_passing_genomes = self_filter_metadata.get(
                    'passing_sample_ids', set())
            other_passing_genomes = other_filter_metadata.get(
                    'passing_sample_ids', set())
            if op_string == '&':
                merged_filter_metadata['passing_sample_ids'] = (
                        self_passing_genomes & other_passing_genomes)
            else:
                merged_filter_metadata['passing_sample_ids'] = (
                        self_passing_genomes | other_passing_genomes)

            # Save the filter metadata.
            new_variant_id_to_metadata_dict[variant.id] = merged_filter_metadata

        return FilterEvalResult(new_variant_set,
                new_variant_id_to_metadata_dict)


FILTER_SCOPE__ALL = 'ALL'
FILTER_SCOPE__ANY = 'ANY'
FILTER_SCOPE__ONLY = 'ONLY'
VALID_FILTER_SCOPES = set([
    FILTER_SCOPE__ALL,
    FILTER_SCOPE__ANY,
    FILTER_SCOPE__ONLY
])

class FilterScope(object):
    """Represents the scope that a filter should be applied over.
    """

    def __init__(self, scope_type, sample_ids):
        """
        Args:
            sample_ids: Set of sample ids.
            scope_type: A scope in VALID_FILTER_SCOPES.
        """
        assert scope_type in VALID_FILTER_SCOPES, "Invalid scope type."

        self.sample_id_set = set(sample_ids)
        self.scope_type = scope_type


    @classmethod
    def parse_sample_ids(clazz, sample_id_string):
        """Turns a comma-separated list of sample uids or names into ids.
        """
        sample_uids_or_names = sample_id_string.split(',')
        sample_uids_or_names = [s.strip() for s in sample_uids_or_names]
        sample_ids = [ExperimentSample.objects.get(uid=uid).id for uid
                in sample_uids_or_names]
        return sample_ids


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

class VariantFilterEvaluator(object):
    """Evaluator for a single scoped expression, e.g. of the form:
        '(position > 5) in ALL(sample1, sample2)'
    """

    def __init__(self, raw_filter_string, ref_genome, scope=None):
        """Constructor.

        Args:
            raw_filter_string: String representing the raw filter.
            ref_genome: ReferenceGenome these variants are relative to.
            scope: Optional FilterScope object which restricts the results
                to the samples according to the semantic setting of the scope.
        """
        # Validation.
        if scope is not None:
            assert isinstance(scope, FilterScope)

        self.raw_filter_string = raw_filter_string
        self.clean_filter_string = raw_filter_string
        self.ref_genome = ref_genome
        self.scope = scope

        # Generator object that provides symbols in alphabetical order.
        self.symbol_maker = symbol_generator()

        # Catch trivial, no filter case.
        if self.clean_filter_string == '':
            self.sympy_representation = ''
        else:
            self._create_symbolic_representation()


    def get_scope_type(self):
        """Returns the scope type.
        """
        if self.scope:
            return self.scope.scope_type
        return None


    def get_scope_sample_id_set(self):
        """Returns the set of sample ids that the scope applies to.
        """
        if self.scope:
            return self.scope.sample_id_set
        return None


    def _create_symbolic_representation(self):
        """Creates a symbolic representation of the query to enable, among
        other things, manipulation with sympy so we can get to disjunctive
        normal form (DNF).
        """
        # Find all the expressions and replace them with symbols.
        self.symbol_to_expression_map = {}

        symbolified_string = self.clean_filter_string
        for regex in [SAMPLE_SCOPE_REGEX, EXPRESSION_REGEX]:
            symbolified_string = self._symbolify_string_for_regex(
                    symbolified_string, regex)

        self.sympy_representation = boolalg.to_dnf(symbolified_string)


    def _symbolify_string_for_regex(self, start_string, regex):
        """Iterates through the string tokenized by the regex and replaces
        the matching tokens with symbols.
        """
        tokens = regex.split(start_string)
        symbolified_tokens = []
        for token in tokens:
            if regex.match(token):
                symbol = self.symbol_maker.next()
                self.symbol_to_expression_map[symbol] = token
                symbolified_tokens.append(symbol)
            else:
                symbolified_tokens.append(token)
        symbolified_string = ''.join(symbolified_tokens)
        return symbolified_string


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
            A FilterEvalResult object.
        """
        return self.evaluate_disjunction(self.sympy_representation)


    def evaluate_disjunction(self, disjunction_clause):
        """Evaluate a disjunction (OR) clause.

        Args:
            disjunction_clause: sympy.logic.boolalg.Or clause, or equivalent.
            scope: FilterScope object.

        Returns:
            A FilterEvalResult object.
        """
        # The disjunction may contain a single condition.
        if not isinstance(disjunction_clause, boolalg.Or):
            return self.evaluate_conjunction(disjunction_clause)
        else:
            result = FilterEvalResult(set(), {})
            for conjunction_clause in disjunction_clause.args:
                result |= self.evaluate_conjunction(conjunction_clause)
            return result


    def evaluate_conjunction(self, conjunction_clause):
        """Evaluate a conjunction condition. That is all symbols are ANDed.

        Args:
            disjunction_clause: sympy.logic.boolalg.And clause, or equivalent.
            scope: FilterScope object.

        Returns:
            A FilterEvalResult object.
        """
        # Iterate through the conditions corresponding to the symbols in the
        # clause and either create Q objects out of them, relative to the
        # Variant model, or save them as key-value models to evaluate in memory
        # after the SQL fetch. Order doesn't matter since this is a conjunction
        # clause.

        def _single_symbol_helper(symbol, filter_eval_results, q_part,
                remaining_triples):
            """Helper method for evaluating a single symbol.
            """
            result = self.handle_single_symbol(symbol)
            if isinstance(result, FilterEvalResult):
                filter_eval_results.append(result)
            elif isinstance(result, Q):
                q_part &= result
            else:
                remaining_triples.append(result)
            return (filter_eval_results, q_part, remaining_triples)

        filter_eval_results = []
        q_part = Q()
        remaining_triples = []

        if not conjunction_clause == '':
            # The conjunction may contain a single condition.
            if not conjunction_clause.is_Boolean:
                (filter_eval_results, q_part, remaining_triples) = (
                        _single_symbol_helper(
                                conjunction_clause, filter_eval_results, q_part,
                                remaining_triples))
            else:
                for symbol in conjunction_clause.args:
                    (filter_eval_results, q_part, remaining_triples) = (
                            _single_symbol_helper(
                                    symbol, filter_eval_results, q_part,
                                    remaining_triples))

        ### Combine any filter results so far. These are probably results
        ### of sub-clauses that are scoped expressions.
        partial_result = None
        if len(filter_eval_results) > 0:
            partial_result = filter_eval_results[0]
            for result in filter_eval_results[1:]:
                partial_result &= result

        ### Now handle the Q part.
        variant_list = self.ref_genome.variant_set.filter(q_part)

        # Make this into a FilterEvalResult object.
        # For now, we just set all the genomes as passing for the conditions
        # so far.
        variant_id_to_metadata_dict = {}
        for variant in variant_list:
            variant_id_to_metadata_dict[variant.id] = {
                'passing_sample_ids': (
                        self.get_sample_id_set_for_variant(variant)),
            }
        q_part_result = FilterEvalResult(set(variant_list),
                variant_id_to_metadata_dict)

        if partial_result is not None:
            partial_result &= q_part_result
        else:
            partial_result = q_part_result

        return self.apply_non_sql_triples_to_query_set(partial_result,
                remaining_triples)


    def handle_single_symbol(self, symbol):
        """Returns one of:
            * A FilterEvalResult object if the symbol represents a scoped
                filter condition.
            * A Django Q object if the symbol represents a condition that
            can be evaluated against the SQL database.
            * A triple of delim, key, value if the condition must be evaluated
                in-memory.
        """
        condition_string = self.get_condition_string_for_symbol(symbol)

        scope_match = SAMPLE_SCOPE_REGEX_NAMED.match(condition_string)
        if scope_match:
            condition_string = scope_match.group('condition')
            scope_type = scope_match.group('scope_type')
            samples_string = scope_match.group('samples')
            sample_ids = FilterScope.parse_sample_ids(samples_string)
            evaluator = VariantFilterEvaluator(condition_string,
                    self.ref_genome, FilterScope(scope_type, sample_ids))
            return evaluator.evaluate()

        (delim, key, value) = _get_delim_key_value_triple(condition_string)
        for key_map in ALL_SQL_KEY_MAP_LIST:
            if key in key_map:
                return _get_django_q_object_for_triple((delim, key, value))
        return (delim, key, value)


    def apply_non_sql_triples_to_query_set(self, filter_eval_result,
            remaining_triples):
        """Applies the remaining condition triples to the query set and
        returns the trimmed down list.
        """
        variant_list = filter_eval_result.variant_set

        variant_id_to_metadata_dict = (
                filter_eval_result.variant_id_to_metadata_dict)

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
                            # No need to update passing sample ids.
                            break

                elif key in VARIANT_EVIDENCE_MAP:
                    samples_passing_for_variant = set()
                    _assert_delim_for_key(VARIANT_EVIDENCE_MAP, delim, key)
                    all_variant_evidence_obj_list = (
                            VariantEvidence.objects.filter(
                                    variant_caller_common_data__in=variant.variantcallercommondata_set.all()))
                    for variant_evidence_obj in all_variant_evidence_obj_list:
                        data_dict = variant_evidence_obj.as_dict()
                        if not data_dict['called']:
                            continue
                        passing = _evaluate_condition_in_triple(
                                data_dict, VARIANT_EVIDENCE_MAP, triple)
                        if passing:
                            samples_passing_for_variant.add(
                                    variant_evidence_obj.experiment_sample.id)

                    # Determine whether the passing samples qualify this
                    # variant as passing the filter, accounting for scope if
                    # applicable.
                    if len(samples_passing_for_variant) > 0:
                        # Compare the passing results to the scope.
                        if self.scope is not None:
                            scope_type = self.get_scope_type()
                            scope_sample_id_set = self.get_scope_sample_id_set()
                            if scope_type == FILTER_SCOPE__ALL:
                                # All passing sample ids must be in the
                                # scope set.
                                intersection = (samples_passing_for_variant &
                                        scope_sample_id_set)
                                if (intersection == scope_sample_id_set):
                                    passing_variant_list.append(variant)
                                    variant_id_to_metadata_dict[variant.id][
                                            'passing_sample_ids'] &= (
                                                    samples_passing_for_variant)
                            elif scope_type == FILTER_SCOPE__ANY:
                                # At least one passing sample id must be in
                                # the scope set.
                                if len(samples_passing_for_variant &
                                        scope_sample_id_set) > 0:
                                    passing_variant_list.append(variant)
                                    variant_id_to_metadata_dict[variant.id][
                                            'passing_sample_ids'] &= (
                                                    samples_passing_for_variant)
                            elif scope_type == FILTER_SCOPE__ONLY:
                                # The passing sample id set must be exactly
                                # the scope set.
                                if (samples_passing_for_variant ==
                                        scope_sample_id_set):
                                    passing_variant_list.append(variant)
                                    variant_id_to_metadata_dict[variant.id][
                                            'passing_sample_ids'] &= (
                                                    samples_passing_for_variant)
                            else:
                                raise AssertionError(
                                        "Unknown scope %s" % scope_type)
                        else:
                            passing_variant_list.append(variant)
                            variant_id_to_metadata_dict[variant.id][
                                    'passing_sample_ids'] &= (
                                            samples_passing_for_variant)
                    else:
                        variant_id_to_metadata_dict[variant.id][
                                'passing_sample_ids'] = set()

                else:
                    raise ParseError(key, 'Unrecognized filter key.')

            # Since we are inside of a conjunction, we only need to check
            # the remaining variants on the next iteration.
            variant_list = passing_variant_list

        return FilterEvalResult(set(variant_list), variant_id_to_metadata_dict)


    def get_sample_id_set_for_variant(self, variant):
        """Returns the set of all ExperimentSamples ids for which there exists a
        relationship to the Variant
        """
        return set([ve.experiment_sample.id for ve in
                VariantEvidence.objects.filter(
                        variant_caller_common_data__variant=variant)])


###############################################################################
# Helper methods
###############################################################################

def _get_delim_key_value_triple(raw_string):
    """Attempt to parse a (delim, key, value) triple out of raw_string."""
    # Remove spaces from the string.
    raw_string = raw_string.replace(' ', '')

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


def _get_django_q_object_for_triple(delim_key_value_triple):
    """Returns a Django Q object for querying against the SNP model for
    the given key_string.

    Args:
        delim_key_value_triple: A tuple representing a single condition.

    Returns a django Q object.
    """
    assert len(delim_key_value_triple) == 3
    (delim, key, value) = delim_key_value_triple

    # Special handling for != delim.
    if delim == '!=':
        postfix = ''
        maybe_not_prefix = '~'
    else:
        postfix = DELIM_TO_Q_POSTFIX[delim]
        maybe_not_prefix = ''

    eval_string = (maybe_not_prefix + 'Q(' + key + postfix + '=' + '"' + value +
            '"' + ')')
    return eval(eval_string)


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
