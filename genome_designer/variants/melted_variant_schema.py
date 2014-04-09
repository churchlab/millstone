"""
Module that builds the melted_variant_schema used to construct the
materialized view of the melted variant data.

Build the schema used to build the materialized view.
"""


class SchemaBuilder(object):
    """Builder object for the schema.

    Usage:
        schema_builder = SchemaBuilder()
        schema_builder.add_melted_variant_field(<source_col_name>,
                <joined_table_col_named>, <is_null_in_variant_to_set_label>,
                <user_queryable>)
        schema_builder.add_melted_variant_field(...)
        schema_builder.add_melted_variant_field(...)
        my_schema = schema_builder.get_schema()
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

# User-queryable schema keys.
# These are all caps to match our code logic of capitalizing the keys in a
# a user's query before handling the filter.
# TODO: Convert all uses of these keys to refer to these constants.
MELTED_SCHEMA_KEY__UID = 'UID'
MELTED_SCHEMA_KEY__POSITION = 'POSITION'
MELTED_SCHEMA_KEY__CHROMOSOME = 'CHROMOSOME'
MELTED_SCHEMA_KEY__REF = 'REF'
MELTED_SCHEMA_KEY__ALT = 'ALT'
MELTED_SCHEMA_KEY__ES_UID = 'EXPERIMENT_SAMPLE_UID'
MELTED_SCHEMA_KEY__ES_LABEL = 'EXPERIMENT_SAMPLE_LABEL'
MELTED_SCHEMA_KEY__VS_UID = 'VARIANT_SET_UID'
MELTED_SCHEMA_KEY__VS_LABEL = 'VARIANT_SET_LABEL'
MELTED_SCHEMA_KEY__VA_ID = 'VA_ID'
MELTED_SCHEMA_KEY__ES_ID = 'ES_ID'
MELTED_SCHEMA_KEY__VE_ID = 'VE_ID'
MELTED_SCHEMA_KEY__VCCD_ID = 'VCCD_ID'

# Used for aggregate total sample count in Postgres query.
CAST_SCHEMA_KEY__TOTAL_SAMPLE_COUNT = 'SAMPLE_COUNT'

SCHEMA_BUILDER = SchemaBuilder()
# SCHEMA_BUILDER.add_melted_variant_field(<source_col_name>,
#         <joined_table_col_named>, <is_null_in_variant_to_set_label>,
#         <user_queryable>)

# Variant
SCHEMA_BUILDER.add_melted_variant_field('main_variant.id', 'id', False, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variant.uid', MELTED_SCHEMA_KEY__UID, False, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variant.position', MELTED_SCHEMA_KEY__POSITION, False, True,
        {'type': 'Integer', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variant.chromosome', MELTED_SCHEMA_KEY__CHROMOSOME, False, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variant.ref_value', MELTED_SCHEMA_KEY__REF, False, True,
        {'type': 'String', 'num': 1})

# VariantAlternate
SCHEMA_BUILDER.add_melted_variant_field('main_variantalternate.id', MELTED_SCHEMA_KEY__VA_ID, False, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variantalternate.alt_value', MELTED_SCHEMA_KEY__ALT, False, True,
        {'type': 'String', 'num': 1})

# Key-value data.
SCHEMA_BUILDER.add_melted_variant_field('main_variantcallercommondata.id', MELTED_SCHEMA_KEY__VCCD_ID, True, False)
SCHEMA_BUILDER.add_melted_variant_field('main_variantevidence.id', MELTED_SCHEMA_KEY__VE_ID, True, False)

# ExperimentSample
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.id', MELTED_SCHEMA_KEY__ES_ID, True, False)
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.uid', MELTED_SCHEMA_KEY__ES_UID, True, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_experimentsample.label', MELTED_SCHEMA_KEY__ES_LABEL, True, True,
        {'type': 'String', 'num': 1})

# VariantSet
SCHEMA_BUILDER.add_melted_variant_field('main_variantset.uid', MELTED_SCHEMA_KEY__VS_UID, False, True,
        {'type': 'String', 'num': 1})
SCHEMA_BUILDER.add_melted_variant_field('main_variantset.label', MELTED_SCHEMA_KEY__VS_LABEL, False, True,
        {'type': 'String', 'num': 1})

# Build the schema.
MELTED_VARIANT_SCHEMA = SCHEMA_BUILDER.get_schema()

# Generate the SELECT clause for building the table.
MATERIALIZED_TABLE_SELECT_CLAUSE_COMPONENTS = []
for schema_obj in MELTED_VARIANT_SCHEMA:
    if (schema_obj['source_col_name'] in [
            'main_variantset.uid',
            'main_variantset.label']):
        MATERIALIZED_TABLE_SELECT_CLAUSE_COMPONENTS.append(
                'array_agg(' + schema_obj['source_col_name'] + ') AS ' +
                schema_obj['joined_table_col_name'])
    else:
        MATERIALIZED_TABLE_SELECT_CLAUSE_COMPONENTS.append(
                schema_obj['source_col_name'] + ' AS ' +
                schema_obj['joined_table_col_name'])
MATERIALIZED_TABLE_SELECT_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_SELECT_CLAUSE_COMPONENTS)

# A GROUP BY, for dealing with repeated variant sets.
MATERIALIZED_TABLE_GROUP_BY_CLAUSE_COMPONENTS = []
for schema_obj in MELTED_VARIANT_SCHEMA:
    if (schema_obj['source_col_name'] in [
            'main_variantset.uid',
            'main_variantset.label']):
        continue
    MATERIALIZED_TABLE_GROUP_BY_CLAUSE_COMPONENTS.append(
            schema_obj['source_col_name'])
MATERIALIZED_TABLE_GROUP_BY_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_GROUP_BY_CLAUSE_COMPONENTS)

# Generate the SELECT clause for the Variant to VariantSet.label view.
# We perform a UNION with this table to ensure that we yield Variants that
# are in a VariantSet without an association with any ExperimentSample.
MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS = []
for schema_obj in MELTED_VARIANT_SCHEMA:
    if schema_obj['is_null_in_variant_to_set_label']:
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS.append(
                'NULL' + ' AS ' + schema_obj['joined_table_col_name'])
    elif (schema_obj['source_col_name'] in [
            'main_variantset.uid',
            'main_variantset.label']):
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS.append(
                'array_agg(' + schema_obj['source_col_name'] + ') AS ' +
                schema_obj['joined_table_col_name'])
    else:
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS.append(
                schema_obj['source_col_name'] + ' AS ' +
                schema_obj['joined_table_col_name'])
MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_VTVS_SELECT_CLAUSE_COMPONENTS)

# A GROUP BY, for dealing with repeated variant sets.
MATERIALIZED_TABLE_VTVS_GROUP_BY_CLAUSE_COMPONENTS = []
for schema_obj in MELTED_VARIANT_SCHEMA:
    if schema_obj['is_null_in_variant_to_set_label']:
        continue
    if (schema_obj['source_col_name'] in [
            'main_variantset.uid',
            'main_variantset.label']):
        continue
    MATERIALIZED_TABLE_VTVS_GROUP_BY_CLAUSE_COMPONENTS.append(
            schema_obj['source_col_name'])
MATERIALIZED_TABLE_VTVS_GROUP_BY_CLAUSE = ', '.join(
        MATERIALIZED_TABLE_VTVS_GROUP_BY_CLAUSE_COMPONENTS)

# Generate the SELECT clause for querying the table.
MATERIALIZED_TABLE_QUERY_SELECT_CLAUSE_COMPONENTS = [
        schema_obj['joined_table_col_name']
        for schema_obj in MELTED_VARIANT_SCHEMA
] + [
    'va_data', # Needed to hard-code showing INFO_EFF_GENE,
    'es_data',
    've_data'
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
