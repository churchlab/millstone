"""
Module that builds the melted_variant_schema used to construct the
materialized view of the melted variant data.

Build the schema used to build the materialized view.
"""


class SchemaBuilder(object):
    """Builder object for the schema.
    """

    def __init__(self):
        self.schema = []

        # Check for duplicates while building.
        self.joined_table_col_name_set = set()

    def add_melted_variant_field(self, source_col_name, joined_table_col_name,
            is_null_in_variant_to_set_label, is_user_queryable,
            query_schema=None):
        assert joined_table_col_name not in self.joined_table_col_name_set
        self.schema.append({
            'source_col_name': source_col_name,
            'joined_table_col_name': joined_table_col_name,
            'is_null_in_variant_to_set_label': is_null_in_variant_to_set_label,
            'is_user_queryable': is_user_queryable,
            'query_schema': query_schema,
        })
        self.joined_table_col_name_set.add(joined_table_col_name)

    def get_schema(self):
        return self.schema


SCHEMA_BUILDER = SchemaBuilder()
# SCHEMA_BUILDER.add_melted_variant_field(<source_col_name>,
#         <joined_table_col_named>, <is_null_in_variant_to_set_label>,
#         <user_queryable>)

# Variant
SCHEMA_BUILDER.add_melted_variant_field('main_variant.id', 'id', False, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.uid', 'uid', False, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variant.position', 'position', False, True,
        {'type': 'Integer', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variant.chromosome', 'chromosome', False, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variant.ref_value', 'ref', False, True,
        {'type': 'String', 'num': 1})

# VariantAlternate
SCHEMA_BUILDER.add_melted_variant_field('main_variantalternate.id', 'va_id', False, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variantalternate.alt_value', 'alt', False, True,
        {'type': 'String', 'num': 1})

# Key-value data.
SCHEMA_BUILDER.add_melted_variant_field('main_variantcallercommondata.id', 'vccd_id', True, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variantevidence.id', 've_id', True, False)

# ExperimentSample
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.id', 'experiment_sample_id', True, False)
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.uid', 'experiment_sample_uid', True, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.label', 'experiment_sample_label', True, True,
        {'type': 'String', 'num': 1})

# VariantSet
SCHEMA_BUILDER.add_melted_variant_field('main_variantset.uid', 'variant_set_uid', False, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variantset.label', 'variant_set_label', False, True,
        {'type': 'String', 'num': 1})

# Build the schema.
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
] + [
    'va_data' # Needed to hard-code showing INFO_EFF_GENE
]

# Map from queryable fields to schema info (e.g. type, num).
MATERIALIZED_TABLE_QUERYABLE_FIELDS_MAP = dict([
        (schema_obj['joined_table_col_name'], schema_obj['query_schema'])
        for schema_obj in MELTED_VARIANT_SCHEMA
        if schema_obj['is_user_queryable']])
# Assert that all user queryable fields have schema defined.
for key, query_schema in MATERIALIZED_TABLE_QUERYABLE_FIELDS_MAP.iteritems():
    assert query_schema is not None, (
            "Missing query schema for queryable %s" % key)
