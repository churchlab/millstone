"""
Methods for filtering over Variants.

Adjusted for the new materialized view strategy. This should end up being less
code than the original implementation; thus it makes sense to re-implement
from scratch.
"""

from django.db import connection
from sympy.logic import boolalg

from variants.common import EXPRESSION_REGEX
from variants.common import SAMPLE_SCOPE_REGEX
from variants.common import GENE_REGEX
from variants.common import SET_REGEX
from variants.common import convert_delim_key_value_triple_to_expr
from variants.common import generate_key_to_materialized_view_parent_col
from variants.common import get_all_key_map
from variants.common import get_delim_key_value_triple
from variants.common import SymbolGenerator
from variants.materialized_view_manager import MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager
from variants.filter_scope import FilterScope


class VariantFilterEvaluator(object):
    """Evaluator for a single scoped expression, e.g. of the form:
        '(position > 5) in ALL(sample1, sample2)'
    """

    def __init__(self, query_args, ref_genome, scope=None):
        """Constructor.

        Args:
            query_args: a dictionary with the following arguments:
                filter_string: String representing the raw filter.
                is_melted: True if melted view, false if cast view.
                sort_by_column: A column name to sort by, or an empty string otherwise.
                count_only: True if only want to return a count
                pagination_start: Offset of the returned query
                pagination_len: Maximum number of returned variants, or -1 for no limit
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

        # Load args into arguments directly
        self.filter_string = query_args.get('filter_string', None)
        self.is_melted = query_args.get('is_melted', True)
        self.sort_by_column = query_args.get('sort_by_column', None)
        self.count_only = query_args.get('count_only', False)
        self.pagination_start = query_args.get('pagination_start', 0)
        self.pagination_len = query_args.get('pagination_len', -1)
        self.visible_key_names = query_args.get('visible_key_names', [])
        self.ref_genome = ref_genome
        self.all_key_map = get_all_key_map(self.ref_genome)
        self.scope = scope

        # Generator object that provides symbols in alphabetical order.
        self.symbol_maker = SymbolGenerator()

        # Catch trivial, no filter case.
        if self.filter_string == '':
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

        symbolified_string = self.filter_string
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
        """Evaluates the given database query.

        Returns:
            A FilterEvalResult object.
        """
        select_clause = self._select_clause()

        cursor = connection.cursor()

        # Minimal sql_statement has select clause.
        sql_statement = 'SELECT %s FROM %s ' % (select_clause,
                self.materialized_view_manager.get_table_name())

        # Maybe add WHERE clause.
        if self.sympy_representation:
            where_clause, where_clause_args = self._where_clause()
            sql_statement += 'WHERE (' + where_clause + ') '
        else:
            where_clause_args = []

        # If cast, need to group by position for array_agg to work.
        if not self.is_melted:
            sql_statement += 'GROUP BY position '

        # Add optional sort clause.
        if self.sort_by_column:
            sql_statement += 'ORDER BY %s ' % self.sort_by_column

        # Add limit and offset clause.
        if self.pagination_len != -1:
            sql_statement += 'LIMIT %d ' % self.pagination_len
        sql_statement += 'OFFSET %d ' % self.pagination_start

        if self.count_only:
            sql_statement = (
                    'SELECT count(*) FROM (' + sql_statement + ') subresult')

        # Execute the query and store the results in hashable representation
        # so that they can be combined through boolean operators with other
        # evaluations.
        cursor.execute(sql_statement, where_clause_args)
        result_list = [dict(zip([col[0] for col in cursor.description], row))
                for row in cursor.fetchall()]
        return result_list

    def _select_clause(self):
        """Determines the SELECT clause for the materialized view.

        If melted, normal select clause with all columns.
        If cast, then array_agg alt and arbitrarily choose one of all
        other fields except for position (by arbitrarily we choose min)

        Returns:
            String representing select clause.
        """
        cols = (MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS +
                self._identify_catch_all_data_fields_to_select())
        if self.is_melted:
            return ', '.join(cols)
        else:
            def fix(column):
                """Helper to fix the aggregate query.

                We GROUP BY position, so all other fields need to be explicitly
                aggregated.
                """
                if column == 'position':
                    return column
                elif column in ['alt', 'va_data', 'vccd_data', 've_data']:
                    return 'array_agg(' + column + ') as ' + column
                else:
                    return 'min(' + column + ') as ' + column
            return ', '.join(map(fix, cols))

    def _identify_catch_all_data_fields_to_select(self):
        """Returns the list of cols to fetch.
        """
        # First build the key to parent map.
        key_to_parent_map = generate_key_to_materialized_view_parent_col(
                self.ref_genome)

        # Now figure out which of the cols to select.
        cols_to_fetch = set()
        for extra_key in self.visible_key_names:
            col = key_to_parent_map.get(extra_key, None)
            if col is not None:
                cols_to_fetch.add(col)
        return list(cols_to_fetch)

    def _where_clause(self):
        # Returns None if no where clause, or
        #   a tuple (where clause, arguments)
        if not self.sympy_representation:
            return None

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
        def _conjuntion_clause(conjunction_clause):
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

            where_clause_conjunctive_expr_list = []
            where_clause_args = []
            for triple in remaining_triples:
                (expr, arg) = convert_delim_key_value_triple_to_expr(triple)
                where_clause_conjunctive_expr_list.append(expr)
                where_clause_args.append(arg)
            return (' AND '.join(where_clause_conjunctive_expr_list),
                    where_clause_args)

        # If only one conjunction clause, evaluate directly
        if not isinstance(self.sympy_representation, boolalg.Or):
            return _conjuntion_clause(self.sympy_representation)
        else:
            conjunction_clause_info = [_conjuntion_clause(conjunction_clause)
                    for conjunction_clause in self.sympy_representation.args]
            where_clause = ' OR '.join(['(' + info[0] + ')'
                for info in conjunction_clause_info])
            args = []
            for info in conjunction_clause_info:
                args.extend(info[1])
            return (where_clause, args)


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

def get_variants_that_pass_filter(query_args, ref_genome):
    """Takes a complete filter string and returns the variants that pass the
    filter.

    We have a hybrid implementation for field values, where some of the field
    values are actually columns in the respective database tables, while others
    are serialized as part of a single key-value string. We provide
    filtering using SQL as much as possible, resorting to per-key filtering
    for special values.

    Args:
        query_args: a dictionary with the following arguments:
            filter_string: String representing the raw filter.
            is_melted: True if melted view, false if cast view.
            sort_by_column: A column name to sort by, or an empty string otherwise.
            count_only: True if only want to return a count
            pagination_start: Offset of the returned query
            pagination_len: Maximum number of returned variants
        ref_genome: The ReferenceGenome this is limited to. This is a hard-coded
            parameter of sorts, since at the least we always limit comparisons
            among Variant objects to those that share a ReferenceGenome.

    Returns:
        List of dictionary objects representing melted Variants.
        See materialized_view_manager.py.
    """
    evaluator = VariantFilterEvaluator(query_args, ref_genome)
    return evaluator.evaluate()
