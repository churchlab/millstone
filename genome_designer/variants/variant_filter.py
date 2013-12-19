"""
Methods for filtering over Variants.

The fields that we allow querying against vary across different models,
and for each model, the field may either be a column in the SQL database, a key
in the catch-all serialized key-value dictionary.

A naive implementation would pull all Variant objects from the SQL database
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
from main.models import Variant
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
from variants.common import Q_OBJECT_TYPE__GLOBAL
from variants.common import Q_OBJECT_TYPE__PER_SAMPLE
from variants.common import VARIANT_ALTERNATE_SQL_KEY_MAP
from variants.common import VARIANT_TABLE_KEY__ID
from variants.common import VARIANT_TABLE_KEY__SAMPLE
from variants.common import assert_delim_for_key
from variants.common import create_initial_filter_eval_result_object
from variants.common import eval_variant_set_filter_expr
from variants.common import evaluate_condition_in_triple
from variants.common import get_all_key_map
from variants.common import get_delim_key_value_triple
from variants.common import get_django_q_object_for_gene_restrict
from variants.common import get_django_q_object_for_triple
from variants.common import get_sample_id_set_for_variant
from variants.common import get_variant_table_column_for_sql_key
from variants.common import SqlReadySymbol
from variants.common import SymbolGenerator
from variants.common import ParseError
from variants.filter_eval_result import FilterEvalResult
from variants.filter_eval_result import metadata_default_dict_factory_fn
from variants.filter_scope import FilterScope


###############################################################################
# Helper objects (also see variants/common.py).
###############################################################################

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
        self.symbol_maker = SymbolGenerator()

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

        This is where the meat of the filtering computation happens.

        Args:
            disjunction_clause: sympy.logic.boolalg.And clause, or equivalent.
            scope: FilterScope object.

        Returns:
            A FilterEvalResult object.
        """
        # On a high level, the algorithm breaks down into doing as much
        # filtering on the SQL side, and then doing the rest in application
        # memory. Fields that we currently handle in memory include those that
        # are part of the catch-all serialized data field, as well as SQL
        # fields which we must then use to establish the specific relationship
        # between Variants and ExperimentSamples.

        # First we bin the individual symbols in the clause according to how
        # they should be handled. The binning decision works according to:
        #     * Symbols that can be handled in SQL are picked out.
        #         * They are all turned into Q objects.
        #         * Those that are relevant to per-sample filtering are saved
        #           for further analysis along with remaining triples below.
        #     * Symbols that are a scoped expression are handled recursively
        #       and the FilterEvalResult object returned.
        #     * Non-SQL triples (part of serialized string) are saved for
        #       for processing in-memory at the end.
        #
        # Order doesn't matter since this is a conjunction (AND) clause.

        # Store the results of binning in these data structures.
        filter_eval_results = []
        sql_ready_symbol_list = []
        remaining_triples = []

        # Perform binning.
        if not conjunction_clause == '':
            # Handle the case where a conjunction might be a single condition.
            if not conjunction_clause.is_Boolean:
                symbol_list = [conjunction_clause]
            else:
                symbol_list = conjunction_clause.args
            for symbol in symbol_list:
                updated_structs = self._single_symbol_mux(
                        symbol, filter_eval_results, sql_ready_symbol_list,
                        remaining_triples)
                (filter_eval_results, sql_ready_symbol_list,
                        remaining_triples) = updated_structs

        # Combine any FilterEvalResults obtained through recursive creation
        # and evaluation of evaluators. These are typically results
        # of sub-clauses that are sample-scoped expressions.
        # TODO: We are also evaluating per-sample SQL queries above, while
        # we should really wait until we apply the global SQL filters, so that
        # it may be possible to search over a smaller Variant space.
        if len(filter_eval_results) > 0:
            partial_result = reduce(lambda accum, iter_val: accum & iter_val,
                    filter_eval_results)
        else:
            partial_result = None

        # Handle global SQL keys. Perform the initial SQL query to constrain
        # the set of Variants that possibly satisfy the query.
        # NOTE: We chain filter calls rather than AND Q objects intentionally.
        #     It is not equivalent to do either.
        global_q_obj_list = [el.q_obj for el in sql_ready_symbol_list
                    if el.semantic_type == Q_OBJECT_TYPE__GLOBAL]
        if len(global_q_obj_list) > 0:
            variant_list = self.ref_genome.variant_set.filter(
                    global_q_obj_list[0])
            for q_obj in global_q_obj_list[1:]:
                variant_list = variant_list.filter(q_obj)
        else:
            variant_list = self.ref_genome.variant_set.all()

        # Create a FilterEvalResult object containing all variants so far
        # and all sample associations.
        q_part_result = create_initial_filter_eval_result_object(variant_list)

        # Now we need to properly handle identify how each ExperimentSample is
        # related to each Variant, if applicable.
        # TODO: Generalize. Repeat for relevant fields.
        relevant_sql_ready_symbol_list = [el for el in sql_ready_symbol_list
                if el.semantic_type == Q_OBJECT_TYPE__PER_SAMPLE]
        if len(relevant_sql_ready_symbol_list):
            ### First, create a joined table with all the data we need.

            # At the least we need the Variant ids, and samples.
            columns_to_get = [
                    VARIANT_TABLE_KEY__ID,
                    VARIANT_TABLE_KEY__SAMPLE]

            # Also request the relevant columns depending on what we are
            # filtering over.
            additional_relevant_column_keys = [
                    get_variant_table_column_for_sql_key(
                            el.delim_key_value_triple[1])
                    for el in relevant_sql_ready_symbol_list]
            columns_to_get += additional_relevant_column_keys

            # Pull the table into memory.
            joined_variant_data_list = Variant.objects.filter(
                    id__in=variant_list).values(*columns_to_get)

            ### Now iterate over the relevant SQL objects and figure out
            ### Variant to ExperimentSample relations.

            for sql_ready_obj in relevant_sql_ready_symbol_list:
                updated_variant_id_to_metadata_dict = defaultdict(
                        metadata_default_dict_factory_fn)

                (delim, key, value) =  sql_ready_obj.delim_key_value_triple

                variant_table_col = get_variant_table_column_for_sql_key(key)

                for variant in joined_variant_data_list:
                    condition_str = (
                            'variant[variant_table_col]' + delim + 'value')
                    if eval(condition_str):
                        updated_variant_id_to_metadata_dict[variant[
                                VARIANT_TABLE_KEY__ID]][
                                        'passing_sample_ids'].add(variant[
                                                VARIANT_TABLE_KEY__SAMPLE])

                # Any variants with 0 associated samples should be removed.
                filtered_variant_list = []
                for variant in variant_list:
                    passing_samples_id_list = updated_variant_id_to_metadata_dict[
                            variant.id]['passing_sample_ids']
                    if len(passing_samples_id_list):
                        filtered_variant_list.append(variant)

                q_part_result &= FilterEvalResult(set(filtered_variant_list),
                        updated_variant_id_to_metadata_dict)

        # Combined with the partial result from above.
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
        elif isinstance(result, SqlReadySymbol):
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
            return eval_variant_set_filter_expr(condition_string,
                    self.ref_genome)

        # Next, check if this is a gene-related expression.
        gene_match = GENE_REGEX.match(condition_string)
        if gene_match:
            return get_django_q_object_for_gene_restrict(condition_string,
                    self.ref_genome)

        # Finally, if here, then this should be a basic, delimiter-separated
        # expression.
        (delim, key, value) = get_delim_key_value_triple(condition_string,
                self.all_key_map)

        # If the key is supported for SQL queries, return the corresponding
        # Q object.
        for key_map in ALL_SQL_KEY_MAP_LIST:
            if key in key_map:
                return get_django_q_object_for_triple((delim, key, value))

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
                    assert_delim_for_key(self.variant_caller_common_map,
                            delim, key)
                    all_common_data_obj = (
                            variant.variantcallercommondata_set.all())
                    # TODO: Figure out semantics of having more than one common
                    # data object.
                    for common_data_obj in all_common_data_obj:
                        if (evaluate_condition_in_triple(
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

        assert_delim_for_key(type_map, delim, key)

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

            if (evaluate_condition_in_triple(data_dict, type_map, triple)):
                samples_passing_for_variant.add(
                        variant_evidence_obj.experiment_sample.id)

        return samples_passing_for_variant


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
