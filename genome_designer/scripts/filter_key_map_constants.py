# Hard-coded keys to be used during parsing.
# Differs from *_SQL_* maps in variants.common as those maps are used
# for determining which keys are supported for SQL-lookups.

SNP_CALLER_COMMON_DATA_HARD_CODED = {
    'CHROM': {'type': 'String', 'num': 1},
    'POS': {'type': 'Integer', 'num': 1},
    'REF': {'type': 'String', 'num': 1},
}

SNP_VARIANT_HARD_CODED = {
    'ALT': {'type': 'String', 'num': -1}
}

SNP_EVIDENCE_HARD_CODED = {
    'GT_TYPE': {'type': 'Integer', 'num': 1},
    'IS_HET': {'type': 'Boolean', 'num': 1}
}

EXPERIMENT_SAMPLE_HARD_CODED = {}

MAP_KEY__VARIANT = 'variant_data'

MAP_KEY__ALTERNATE = 'snp_alternate_data'

MAP_KEY__COMMON_DATA = 'snp_caller_common_data'

MAP_KEY__EVIDENCE = 'snp_evidence_data'

MAP_KEY__EXPERIMENT_SAMPLE = 'experiment_sample_data'