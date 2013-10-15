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

from collections import defaultdict
from django.db.models import Q
from sympy.logic import boolalg

from main.models import ExperimentSample
from main.models import Region
from main.models import VariantAlternate
from main.models import VariantEvidence
from variants.common import ALL_SQL_KEY_MAP_LIST
from variants.common import DELIM_TO_Q_POSTFIX
from variants.common import EXPRESSION_REGEX
from variants.common import SAMPLE_SCOPE_REGEX
from variants.common import SAMPLE_SCOPE_REGEX_NAMED
from variants.common import GENE_REGEX
from variants.common import GENE_REGEX_NAMED
from variants.common import SET_REGEX
from variants.common import SET_REGEX_NAMED
from variants.common import TYPE_TO_SUPPORTED_OPERATIONS
from variants.common import VARIANT_ALTERNATE_SQL_KEY_MAP
from variants.common import get_all_key_map
from variants.common import get_delim_key_value_triple
from variants.common import ParseError
from variants.filter_eval_result import FilterEvalResult
from variants.filter_scope import FilterScope


###############################################################################
# Constants
###############################################################################

# Q-object types. This is used to determine whether we need to perform a
# query per sample.
Q_OBJECT_TYPE__GLOBAL = 'global'
Q_OBJECT_TYPE__PER_SAMPLE = 'per_sample'


###############################################################################
# Helper objects.
###############################################################################

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
        self.all_key_map = get_all_key_map(self.ref_genome)
        self.variant_caller_common_map = (
                self.ref_genome.get_variant_caller_common_map())
        self.variant_alternate_map = (
                self.ref_genome.get_variant_alternate_map())
        self.variant_evidence_map = (
                self.ref_genome.get_variant_evidence_map())
        self.scope = scope

        # Generator object that provides symbols in alphabetical order.
        self.symbol_maker = symbol_generator()

        # Catch trivial, no filter case.
        if self.clean_filter_string == '':
            self.sympy_representation = ''
        else:
            self._create_symbolic_representation()


    def _create_symbolic_representation(self):
        """Creates a symbolic representation of the query to enable, among
        other things, manipulation with sympy so we can get to disjunctive
        normal form (DNF).
        """
        # Find all the expressions and replace them with symbols.
        self.symbol_to_expression_map = {}

        symbolified_string = self.clean_filter_string
        for regex in [SAMPLE_SCOPE_REGEX, EXPRESSION_REGEX, SET_REGEX,
                GENE_REGEX]:
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

        filter_eval_results = []
        q_list = []
        remaining_triples = []

        if not conjunction_clause == '':
            # The conjunction may contain a single condition.
            if not conjunction_clause.is_Boolean:
                (filter_eval_results, q_list, remaining_triples) = (
                        self._single_symbol_mux(
                                conjunction_clause, filter_eval_results, q_list,
                                remaining_triples))
            else:
                for symbol in conjunction_clause.args:
                    (filter_eval_results, q_list, remaining_triples) = (
                            self._single_symbol_mux(
                                    symbol, filter_eval_results, q_list,
                                    remaining_triples))

        # Combine any filter results so far. These are probably results
        # of sub-clauses that are scoped expressions.
        partial_result = None
        if len(filter_eval_results) > 0:
            partial_result = filter_eval_results[0]
            for result in filter_eval_results[1:]:
                partial_result &= result

        # TODO: Handle evaluating sets for melted view.

        ### Now handle the Q list.

        # NOTE: Previously, we were applying boolean operators
        #     among Q objects, but the ORM mapping fails in cases where you
        #     want to handle logic among VariantSet membership. Thus we chain
        #     filter() calls instead.
        if len(q_list) > 0:
            actual_q_obj_list = [el[1] for el in q_list]
            variant_list = self.ref_genome.variant_set.filter(
                    actual_q_obj_list[0])
            for q_obj in actual_q_obj_list[1:]:
                variant_list = variant_list.filter(q_obj)
        else:
            variant_list = self.ref_genome.variant_set.all()

        # Make this into a FilterEvalResult object, which requires creating the
        # map from variant id to sample ids passing for that variant.

        # First we set all samples passing.
        variant_id_to_metadata_dict = {}
        for variant in variant_list:
            variant_id_to_metadata_dict[variant.id] = {
                'passing_sample_ids': get_sample_id_set_for_variant(variant),
            }

        # Now, for the Q objects that may give different results per-sample,
        # we need to test each specific sample against the query.
        # NOTE: As a first stab, we are going to only make this work for
        #     alt_values.
        # TODO: Come up with a more general implementation that works for sets
        #     and genes.
        per_alt_q_obj_list = [el[1] for el in q_list
                if el[0] == Q_OBJECT_TYPE__PER_SAMPLE]
        if len(per_alt_q_obj_list) > 0:
            assert len(per_alt_q_obj_list) == 1, "Only support 1 right now."
            q_obj = per_alt_q_obj_list[0]
            for variant in variant_list:
                # HACK: This is super ghetto but just trying to get the test to pass
                # for now.
                per_sample_field_name = q_obj.children[0][0].split('__')[-1]
                combined_q = eval('Q(variantalternate_set__' +
                        per_sample_field_name + '=' + '"' +
                        q_obj.children[0][1] + '"' + ')')
                # For now, we only handle the alt_values key.
                passing_evidence_obj_list = (
                        VariantEvidence.objects.filter(combined_q))
                passing_sample_id_set = set([evidence.experiment_sample.id for
                        evidence in passing_evidence_obj_list])
                variant_id_to_metadata_dict[variant.id] = {
                    'passing_sample_ids': passing_sample_id_set,
                }

        # Now build the filter result object.
        q_part_result = FilterEvalResult(set(variant_list),
                variant_id_to_metadata_dict)

        if partial_result is not None:
            partial_result &= q_part_result
        else:
            partial_result = q_part_result

        return self.apply_non_sql_triples_to_query_set(partial_result,
                remaining_triples)


    def _single_symbol_mux(self, symbol, filter_eval_results, q_list,
            remaining_triples):
        """Helper method for evaluating a single symbol.
        """
        result = self.handle_single_symbol(symbol)
        if isinstance(result, FilterEvalResult):
            filter_eval_results.append(result)
        elif isinstance(result, tuple) and isinstance(result[1], Q):
            q_list.append(result)
        else:
            remaining_triples.append(result)
        return (filter_eval_results, q_list, remaining_triples)


    def handle_single_symbol(self, symbol):
        """Returns one of:
            * A FilterEvalResult object if the symbol represents a scoped
                filter condition.
            * A Tuple pair (Q_OBJ_TYPE, Django Q object), if the symbol
                represents a condition that
            can be evaluated against the SQL database.
            * A triple of delim, key, value if the condition must be evaluated
                in-memory.
        """
        condition_string = self.get_condition_string_for_symbol(symbol)

        ### First we check for complex expressions (e.g. scoped, set, etc.)

        # First, check whether this expression is contained within a scope.
        scope_match = SAMPLE_SCOPE_REGEX_NAMED.match(condition_string)
        if scope_match:
            condition_string = scope_match.group('condition')
            scope_type = scope_match.group('scope_type')
            samples_string = scope_match.group('samples')
            sample_ids = FilterScope.parse_sample_ids(samples_string)
            evaluator = VariantFilterEvaluator(condition_string,
                    self.ref_genome, FilterScope(scope_type, sample_ids))
            return evaluator.evaluate()

        # Next, check if this is a set-related expression.
        set_match = SET_REGEX.match(condition_string)
        if set_match:
            return _get_django_q_object_for_set_restrict(condition_string)

        # Next, check if this is a gene-related expression.
        gene_match = GENE_REGEX.match(condition_string)
        if gene_match:
            return _get_django_q_object_for_gene_restrict(condition_string,
                    self.ref_genome)

        ### If we're here, then the condition should be a basic,
        ### delimiter-separated expression.

        # Finally, if here, then this should be a basic, delimiter-separated
        # expression.
        (delim, key, value) = get_delim_key_value_triple(condition_string,
                self.all_key_map)

        # If the key is supported for SQL queries, return the corresponding
        # Q object.
        for key_map in ALL_SQL_KEY_MAP_LIST:
            if key in key_map:
                return _get_django_q_object_for_triple((delim, key, value))

        # Otherwise just return the triple to be evaluated separately.
        return (delim, key, value)


    def apply_non_sql_triples_to_query_set(self, filter_eval_result,
            remaining_triples):
        """Applies the remaining condition triples to the query set and
        returns the trimmed down list.
        """
        # Parse the given FilterEvalResult object.
        variant_list = filter_eval_result.variant_set
        variant_id_to_metadata_dict = (
                filter_eval_result.variant_id_to_metadata_dict)

        for triple in remaining_triples:
            (delim, key, value) = triple

            passing_variant_list = []
            for variant in variant_list:

                # TODO: Currently we are treating alternate keys ENTIRELY as if
                # they were specific to samples. This fine EXCEPT in cases where
                # there exists an INFO value which none of the samples have. The
                # expected result would be to return the variant but no samples,
                # which will not happen if we do it this way.

                # First make sure this is a valid key.
                if not (key in self.variant_caller_common_map or
                        key in self.variant_alternate_map or 
                        key in self.variant_evidence_map):
                    raise ParseError(key, 'Unrecognized filter key.')

                if key in self.variant_caller_common_map:
                    _assert_delim_for_key(self.variant_caller_common_map,
                            delim, key)
                    all_common_data_obj = (
                            variant.variantcallercommondata_set.all())
                    # TODO: Figure out semantics of having more than one common
                    # data object.
                    for common_data_obj in all_common_data_obj:
                        if (_evaluate_condition_in_triple(
                                common_data_obj.as_dict(),
                                self.variant_caller_common_map,
                                triple)):
                            passing_variant_list.append(variant)
                            # No need to update passing sample ids.
                            break

                else: # (if key is per-sample or per-alternate)
                    samples_passing_for_variant = (
                            self.get_samples_passing_for_evidence_or_alternate(
                                    variant, triple))

                    # Determine whether the passing samples qualify this
                    # variant as passing the filter, accounting for scope if
                    # applicable.
                    if len(samples_passing_for_variant) > 0:
                        # Either there is no scope, or samples must pass
                        # the scope.
                        if (self.scope is None or
                                self.scope.do_passing_samples_satisfy_scope(
                                        samples_passing_for_variant)):
                            passing_variant_list.append(variant)
                            variant_id_to_metadata_dict[variant.id][
                                    'passing_sample_ids'] &= (
                                            samples_passing_for_variant)
                    else:
                        # NOTE: The reason we still store the empty set is
                        #     so that we can perform mergers of result objects
                        #     (I think).
                        variant_id_to_metadata_dict[variant.id][
                            'passing_sample_ids'] = set()

            # Since we are inside of a conjunction, we only need to check
            # the remaining variants on the next iteration.
            variant_list = passing_variant_list

        return FilterEvalResult(set(variant_list), variant_id_to_metadata_dict)


    def get_samples_passing_for_evidence_or_alternate(self, variant, triple):

        (delim, key, value) = triple
        samples_passing_for_variant = set()

        # Use the appropriate type map if this is a per-alt key.
        if key in self.variant_alternate_map:
            type_map = self.variant_alternate_map
        elif key in self.variant_evidence_map:
            type_map = self.variant_evidence_map
        else:
            raise InputError('Key passed is not in evidence or alternate map.')

        _assert_delim_for_key(type_map, delim, key)

        # Check each VariantEvidence object.
        all_variant_evidence_obj_list = (
                VariantEvidence.objects.filter(
                        variant_caller_common_data__variant=variant))

        for variant_evidence_obj in all_variant_evidence_obj_list:
            data_dict = variant_evidence_obj.as_dict()

            if not data_dict['called']:
                continue

            # For VariantAlternate (per-alt) keys, map the sample's alleles
            # onto items in the list. For instance, if a sample has
            # a genotype of 1/1, then it's items will be the first
            # allele in the -1 list of INFO_EFF_* fields.
            if type_map == self.variant_alternate_map:
                alts = VariantAlternate.objects.filter(
                    variant=variant,
                    variantevidence=variant_evidence_obj)

                data_dict = defaultdict(list)
                for alt in alts:
                    alt_dict = alt.as_dict()
                    [data_dict[k].append(alt_dict[k]) for k in alt_dict.keys()]

            if (_evaluate_condition_in_triple(data_dict, type_map, triple)):
                samples_passing_for_variant.add(
                        variant_evidence_obj.experiment_sample.id)

        return samples_passing_for_variant


###############################################################################
# Helper methods
###############################################################################

def get_sample_id_set_for_variant(variant):
    """Returns the set of all ExperimentSamples ids for which there exists a
    relationship to the Variant
    """
    return set([ve.experiment_sample.id for ve in
            VariantEvidence.objects.filter(
                    variant_caller_common_data__variant=variant)])


def _get_django_q_object_for_triple(delim_key_value_triple):
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


def _get_django_q_object_for_set_restrict(set_restrict_string):
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


def _get_django_q_object_for_gene_restrict(gene_restrict_string, ref_genome):
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


def _evaluate_condition_in_triple(data_map, type_map, triple, idx=None):
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
            evaluations.append(_evaluate_condition_in_triple(
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
