"""
Methods for filtering over Variants.

Adjusted for the new materialized view strategy. This should end up being less
code than the original implementation; thus it makes sense to re-implement
from scratch.
"""

from django.db import connection
from sympy.logic import boolalg

from variants.common import ALL_SQL_KEY_MAP_LIST
from variants.common import DELIM_TO_Q_POSTFIX
from variants.common import EXPRESSION_REGEX
from variants.common import SAMPLE_SCOPE_REGEX
from variants.common import SAMPLE_SCOPE_REGEX_NAMED
from variants.common import GENE_REGEX
from variants.common import GENE_REGEX_NAMED
from variants.common import SET_REGEX
from variants.common import SET_REGEX_NAMED
from variants.common import dictfetchall
from variants.common import convert_delim_key_value_triple_to_expr
from variants.common import get_all_key_map
from variants.common import get_delim_key_value_triple
from variants.common import get_django_q_object_for_triple
from variants.common import hashablefetchall
from variants.common import SqlReadySymbol
from variants.common import SymbolGenerator
from variants.materialized_view_manager import MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager
from variants.filter_eval_result import FilterEvalResult
from variants.filter_scope import FilterScope


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
        # Manager for making queries to the materialized view table
        # for this ReferenceGenome.
        self.materialized_view_manager = MeltedVariantMaterializedViewManager(
                ref_genome)
        self.materialized_view_manager.create_if_not_exists()

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

        Returns:
            A FilterEvalResult object.
        """
        # The disjunction may contain a single condition.
        if not isinstance(disjunction_clause, boolalg.Or):
            return self.evaluate_conjunction(disjunction_clause)
        else:
            result = FilterEvalResult(set())
            for conjunction_clause in disjunction_clause.args:
                result |= self.evaluate_conjunction(conjunction_clause)
            return result


    def evaluate_conjunction(self, conjunction_clause):
        """Evaluate a conjunction condition. That is all symbols are ANDed.

        This is where the meat of the filtering computation happens.

        Args:
            conjunction_clause: sympy.logic.boolalg.And clause, or equivalent.

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

        # Bin the parts of the query, potentially making recursive calls
        # if we have scoped filters (e.g. per-sample).
        if not conjunction_clause == '':
            # Handle the case where a conjunction might be a single condition.
            if not conjunction_clause.is_Boolean:
                symbol_list = [conjunction_clause]
            else:
                symbol_list = conjunction_clause.args

            # Bin.
            for symbol in symbol_list:
                updated_structs = self._single_symbol_mux(
                        symbol, filter_eval_results, sql_ready_symbol_list,
                        remaining_triples)
                (filter_eval_results, sql_ready_symbol_list,
                        remaining_triples) = updated_structs

        # Handle global SQL keys. Perform the initial SQL query to constrain
        # the results.
        cursor = connection.cursor()
        sql_statement = (
                'SELECT %s '
                'FROM %s '
                % (MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE,
                        self.materialized_view_manager.get_table_name())
        )

        # Maybe add WHERE clause.
        where_clause_conjunctive_expr_list = []
        where_clause_args = []
        for triple in remaining_triples:
            (expr, arg) = convert_delim_key_value_triple_to_expr(triple)
            where_clause_conjunctive_expr_list.append(expr)
            where_clause_args.append(arg)
        where_clause_content = ' AND '.join(where_clause_conjunctive_expr_list)
        if where_clause_content:
            sql_statement += 'WHERE (' + where_clause_content + ') '

        # Execute the query and store the results in hashable representation
        # so that they can be combined through boolean operators with other
        # evaluations.
        cursor.execute(sql_statement, where_clause_args)
        result_list = hashablefetchall(cursor)
        return FilterEvalResult(result_list)


    def _single_symbol_mux(self, symbol, filter_eval_results, q_list,
            remaining_triples):
        """Helper method for evaluating a single symbol.
        """
        result = self._handle_single_symbol(symbol)
        remaining_triples.append(result)
        return (filter_eval_results, q_list, remaining_triples)


    def _handle_single_symbol(self, symbol):
        """Returns one of:
            * A FilterEvalResult object if the symbol represents a scoped
                filter condition.
            * A Tuple pair (Q_OBJ_TYPE, Django Q object), if the symbol
                represents a condition that
            can be evaluated against the SQL database.
            * A triple of delim, key, value if the condition must be evaluated
                in-memory.
        """
        # Look up the expression that the symbol represents.
        condition_string = self.get_condition_string_for_symbol(symbol)

        # Finally, if here, then this should be a basic, delimiter-separated
        # expression.
        (delim, key, value) = get_delim_key_value_triple(condition_string,
                self.all_key_map)
        return (delim, key, value)


###############################################################################
# Main client method.
###############################################################################

def get_variants_that_pass_filter(filter_string, ref_genome):
    """Takes a complete filter string and returns the variants that pass the
    filter.

    We have a hybrid implementation for field values, where some of the field
    values are actually columns in the respective database tables, while others
    are serialized (pickled) as part of a single key-value string. We provide
    filtering using SQL as much as possible, resorting to per-key filtering
    for special values.

    Args:
        filter_string: Query string from the user.
        ref_genome: The ReferenceGenome this is limited to. This is a hard-coded
            parameter of sorts, since at the least we always limit comparisons
            among Variant objects to those that share a ReferenceGenome.

    Returns:
        List of dictionary objects representing melted Variants.
        See materialized_view_manager.py.
    """
    evaluator = VariantFilterEvaluator(filter_string, ref_genome)
    return evaluator.evaluate()
