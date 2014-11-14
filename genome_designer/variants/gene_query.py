"""
Code for querying by genes.

We leverage the MATERIALIZED VIEW for Variants.
"""

from django.db import connection

from variants.materialized_view_manager import MeltedVariantMaterializedViewManager


def lookup_genes(alignment_group):
	"""Looks up Genes.

	Returns list of dictionaries with keys:
		* gene
		* num_variants
	"""
	# Perform the query against the melted variant view.
	materialized_view_manager = MeltedVariantMaterializedViewManager(
		alignment_group.reference_genome)
	materialized_view_manager.create_if_not_exists_or_invalid()

	# Build up the sql statement in parts.

	# Select the gene data and relevant counts.
	select_clause = (
		"va_data->>'INFO_EFF_GENE' AS gene, "
		"COUNT(DISTINCT position) AS num_variants ")

	# Start building the sql statement.
	sql_statement = 'SELECT %s FROM %s ' % (select_clause,
		materialized_view_manager.get_table_name())

	# Add the where clause.
	where_clause_gene_part = "((va_data->>'INFO_EFF_GENE'::text) IS NOT NULL) "
	where_clause_alignment_group_part = (
			'AG_ID = {ag_id} OR AG_ID IS NULL'.format(
					ag_id=alignment_group.id))
	where_clause = '({ag_part}) AND ({gene_part})'.format(
			ag_part=where_clause_alignment_group_part,
			gene_part=where_clause_gene_part)
	sql_statement += 'WHERE (' + where_clause + ') '

	# Finally group by gene.
	sql_statement += 'GROUP BY gene '

	# Execute the query.
	cursor = connection.cursor()
	cursor.execute(sql_statement)

    # Column header data.
	col_descriptions = [col[0].upper() for col in cursor.description]

	return [dict(zip(col_descriptions, row)) for row in cursor.fetchall()]
