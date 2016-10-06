REF_GENOME_UID = '92c717c3'

cols = [
    'id',
    'uid',
    'position',
    'chromosome',
    'ref',
    'va_id',
    'alt',
    'ag_id',
    'vccd_id',
    've_id',
    'es_id',
    'experiment_sample_uid',
    'experiment_sample_label',
    'variant_set_uid',
    'variant_set_label',
    'va_data',
    'vccd_data',
    've_data',
    'es_data',
]

with open('millstone_create_indeces_cmd.txt', 'w') as fh:
    for c in cols:
        fh.write(
            'CREATE INDEX mmv_{ref_uid}_{col} on materialized_melted_variant_{ref_uid} ({col});\n'.format(
                ref_uid=REF_GENOME_UID,
                col=c))
