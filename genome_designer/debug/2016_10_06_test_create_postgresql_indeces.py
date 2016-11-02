"""Script that generates Postgres commands to create indeces on materialized view.

TODO(gleb): Automate this functionality.
"""

import sys


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


json_keys = {
        'va_data': [
                'INFO_AF',
                'INFO_AO',
                'INFO_DPRA',
                'INFO_EFF_CLASS',
                'INFO_EFF_CODING',
                'INFO_EFF_EFFECT',
                'INFO_EFF_GENE',
                'INFO_EFF_IMPACT',
                'INFO_LEN',
                'INFO_SV_LEN',
                'INFO_TYPE'
        ],

        'vccd_data': [
            'INFO_DP',
            'INFO_METHOD',
            'INFO_NS',
            'INFO_SV_TYPE',
            'IS_SV',
        ],

        've_data': [
                'AF',
                'DP',
                'GT_TYPE',
                'IS_HET',
                'RO',
                'AO',
        ],

        'es_data': [
                'SAMPLE_FITNESS',
                'SAMPLE_LINE',
                'SAMPLE_PARENTS'
        ],
}


if len(sys.argv) < 2:
    print 'Usage: python 2016_10_06_test_create_postgresql_indeces.py <ref_genome_uid>'
    sys.exit(0)
else:
    ref_genome_uid_list = sys.argv[1:]
    print 'Preparing for ref genome uids %s' % ref_genome_uid_list

with open('millstone_create_indeces_cmd.txt', 'w') as fh:
    for ref_genome_uid in ref_genome_uid_list:
         for c in cols:
             fh.write(
                 'CREATE INDEX mmv_{ref_uid}_{col} on materialized_melted_variant_{ref_uid} ({col});\n'.format(
                     ref_uid=ref_genome_uid,
                     col=c))
         for data_key, subkey_list in json_keys.iteritems():
             for subkey in subkey_list:
                fh.write(
                    'CREATE INDEX mmv_{ref_uid}_{data_key}_{subkey} on materialized_melted_variant_{ref_uid}(({data_key}->>\'{subkey}\'));\n'.format(
                        ref_uid=ref_genome_uid,
                        data_key=data_key,
                        subkey=subkey))

print 'Done.'
