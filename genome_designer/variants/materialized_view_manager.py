"""
Manages the Materialized view of the Variant data for filtering.
"""

from django.db import connection
from django.db import transaction

class AbstractMaterializedViewManager(object):
    """Base class for object acting as wrapper for a Postgresql materialized
    view.
    """

    def __init(self):
        self.view_table_name = None


    def create(self):
        """Creates the materialized view in the Postgresql DB.
        """
        # Drop the existing table.
        self.drop()

        # Delegate to child class.
        self.create_internal()


    def create_internal(self):
        """Creates the materialized view in the Postgresql DB.

        Child classes should implement.
        """
        pass


    def refresh(self):
        """Refreshes the view.
        """
        refresh_statement = 'REFRESH MATERIALIZED VIEW %s ' % (
                self.view_table_name)
        self.cursor.execute(refresh_statement)


    def drop(self):
        """Drops the materialized view in the Postgresql DB.
        """
        assert self.view_table_name
        drop_sql_statement = "DROP MATERIALIZED VIEW IF EXISTS %s" % (
                self.view_table_name,)
        self.cursor.execute(drop_sql_statement)
        transaction.commit_unless_managed()


# Build the schema used to build the materialized view.
# All of this happens on the first module import.
class SchemaBuilder(object):
    def __init__(self):
        self.schema = []

        # Check for duplicates while building.
        self.joined_table_col_name_set = set()

    def add_melted_variant_field(self, source_col_name, joined_table_col_name,
            is_null_in_variant_to_set_label, is_user_queryable):
        assert joined_table_col_name not in self.joined_table_col_name_set
        self.schema.append({
            'source_col_name': source_col_name,
            'joined_table_col_name': joined_table_col_name,
            'is_null_in_variant_to_set_label': is_null_in_variant_to_set_label,
            'is_user_queryable': is_user_queryable,
        })
        self.joined_table_col_name_set.add(joined_table_col_name)

    def get_schema(self):
        return self.schema

SCHEMA_BUILDER = SchemaBuilder()
# SCHEMA_BUILDER.add_melted_variant_field(<source_col_name>,
#         <joined_table_col_named>, <is_null_in_variant_to_set_label>,
#         <user_queryable>)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.id', 'id', False, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.uid', 'uid', False, True)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.position', 'position', False, True)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.chromosome', 'chromosome', False, True)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.ref_value', 'ref', False, True)
SCHEMA_BUILDER.add_melted_variant_field('main_variantalternate.alt_value', 'alt', False, True)
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.id', 'experiment_sample_id', True, False)
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.uid', 'experiment_sample_uid', True, True)
SCHEMA_BUILDER.add_melted_variant_field('main_variantset.uid', 'variant_set_uid', False, True)
SCHEMA_BUILDER.add_melted_variant_field('main_variantset.label', 'variant_set_label', False, True)
MELTED_VARIANT_SCHEMA = SCHEMA_BUILDER.get_schema()

# Generate the SELECT clause for building the table.
MATERIALIZED_TABLE_SELECT_CLAUSE_COMPONENTS = [
        schema_obj['source_col_name'] + ' AS ' + schema_obj['joined_table_col_name']
        for schema_obj in MELTED_VARIANT_SCHEMA
]
MATERIALIZED_TABLE_SELECT_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_SELECT_CLAUSE_COMPONENTS)

# Generate the SELECT clause for the Variant to VariantSet.label view.
# We perform a UNION with this table to ensure that we yield Variants that
# are in a VariantSet without an association with any ExperimentSample.
MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS = []
for schema_obj in MELTED_VARIANT_SCHEMA:
    if schema_obj['is_null_in_variant_to_set_label']:
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS.append(
                'NULL' + ' AS ' + schema_obj['joined_table_col_name'])
    else:
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS.append(
                schema_obj['source_col_name'] + ' AS ' +
                schema_obj['joined_table_col_name'])
MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS)

# Generate the SELECT clause for querying the table.
MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS = [
        schema_obj['joined_table_col_name']
        for schema_obj in MELTED_VARIANT_SCHEMA
]
MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS)


class MeltedVariantMaterializedViewManager(AbstractMaterializedViewManager):
    """Interface for objects providing a wrapper for a Postgresql materialized
    view.
    """

    def __init__(self, reference_genome):
        self.reference_genome = reference_genome
        self.view_table_name = self.get_table_name()
        self.cursor = connection.cursor()


    def get_table_name(self):
        """Get the name of the underlying SQL table.
        """
        return 'materialized_melted_variant_' + self.reference_genome.uid


    def check_table_exists(self):
        """Check if the table exists.

        NOTE: Figured out the raw sql query by running psql with -E flag
        and then calling \d. The -E flag causes the raw sql of the commands
        to be shown.
        """
        raw_sql = (
            'SELECT c.relname '
            'FROM pg_catalog.pg_class c '
            'WHERE c.relkind=%s AND c.relname=%s '
        )
        self.cursor.execute(raw_sql, ('m', self.view_table_name))
        return bool(self.cursor.fetchone())


    def create_internal(self):
        """Override.
        """
        create_sql_statement = (
            'CREATE MATERIALIZED VIEW %s AS '
                '(SELECT %s FROM main_variant '
                    'INNER JOIN main_variantcallercommondata ON (main_variant.id = main_variantcallercommondata.variant_id) '
                    'INNER JOIN main_variantevidence ON (main_variantcallercommondata.id = main_variantevidence.variant_caller_common_data_id) '
                    'INNER JOIN main_experimentsample ON (main_variantevidence.experiment_sample_id = main_experimentsample.id) '

                    # VariantSet
                    'LEFT JOIN main_varianttovariantset_sample_variant_set_association ON ('
                            'main_experimentsample.id = main_varianttovariantset_sample_variant_set_association.experimentsample_id) '
                    'LEFT JOIN main_varianttovariantset ON ('
                            'main_varianttovariantset_sample_variant_set_association.varianttovariantset_id = main_varianttovariantset.id AND '
                            'main_varianttovariantset.variant_id = main_variant.id) '
                    'LEFT JOIN main_variantset ON main_varianttovariantset.variant_set_id = main_variantset.id '

                    # VariantAlternate
                    'LEFT JOIN main_variantevidence_variantalternate_set ON ('
                            'main_variantevidence.id = main_variantevidence_variantalternate_set.variantevidence_id) '
                    'LEFT JOIN main_variantalternate ON main_variantevidence_variantalternate_set.variantalternate_id = main_variantalternate.id '
                'WHERE (main_variant.reference_genome_id = %d)) '
                'UNION '
                '(SELECT %s FROM main_variant '
                    'INNER JOIN main_variantalternate ON main_variantalternate.variant_id = main_variant.id '
                    'INNER JOIN main_varianttovariantset ON main_variant.id = main_varianttovariantset.variant_id '
                    'INNER JOIN main_variantset ON main_varianttovariantset.variant_set_id = main_variantset.id '
                'WHERE (main_variant.reference_genome_id = %d)) '
                'ORDER BY position, experiment_sample_uid DESC '
            % (self.view_table_name, MATERIALIZED_TABLE_SELECT_CLAUSE, self.reference_genome.id,
                    MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE, self.reference_genome.id)
        )
        self.cursor.execute(create_sql_statement)
        transaction.commit_unless_managed()

    def create_if_not_exists(self):
        if not self.check_table_exists():
            self.create()
