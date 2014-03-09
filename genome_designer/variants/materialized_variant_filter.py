"""
Methods for filtering over Variants.

Adjusted for the new materialized view strategy. This should end up being less
code than the original implementation; thus it makes sense to re-implement
from scratch.
"""

from django.db import connection

from variants.materialized_view_manager import MATERIALIZED_TABLE_MINIMAL_QUERY_SELECT_CLAUSE
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
        self.filter_string = query_args['filter_string']
        self.is_melted = query_args['is_melted']
        self.sort_by_column = query_args['sort_by_column']
        self.count_only = query_args['count_only']
        self.pagination_start = query_args['pagination_start']
        self.pagination_len = query_args['pagination_len']
        self.ref_genome = ref_genome
        self.scope = scope


    def evaluate(self):
        """Evaluates the given database query.

        Returns:
            A FilterEvalResult object.
        """
        # Create select clause. If melted, normal select clause with all columns.
        # If cast, then array_agg alt and arbitrarily choose one of all
        # other fields except for position (by arbitrarily we choose min)
        if self.is_melted:
            select_clause = MATERIALIZED_TABLE_MINIMAL_QUERY_SELECT_CLAUSE
        else:
            columns = MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS
            def fix(column):
                if column == 'position':
                    return column
                elif column == 'alt':
                    return 'array_agg(alt) as alt'
                else:
                    return 'min(' + column + ') as ' + column
            select_clause = ', '.join(map(fix, columns))

        # Get table name
        table_name = self.materialized_view_manager.get_table_name()

        # Handle global SQL keys. Perform the initial SQL query to constrain
        # the results.
        cursor = connection.cursor()
        sql_statement = 'SELECT %s FROM %s ' % (select_clause, table_name)

        # Maybe add WHERE clause.
        if self.filter_string:
            sql_statement += 'WHERE (' + self.filter_string + ') '

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
            sql_statement = 'SELECT count(*) FROM (' + sql_statement + ') subresult';

        # Execute the query and store the results in hashable representation
        # so that they can be combined through boolean operators with other
        # evaluations.
        cursor.execute(sql_statement)
        result_list = [dict(zip([col[0] for col in cursor.description], row))
                for row in cursor.fetchall()]
        return result_list


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
