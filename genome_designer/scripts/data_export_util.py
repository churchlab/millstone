"""
Methods for exporting data.
"""

import csv

from django.db import connection

from variants.common import dictfetchall
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager


def export_melted_variant_view(ref_genome, variant_id_list, csvfile):
    """Exports Variants as a flat .csv file.

    Args:
        variant_id_list: List of variant ids to query. NOTE: Unused.
        csfile: File object to which we write the result.
    """
    mvm = MeltedVariantMaterializedViewManager(ref_genome)

    # Perform the query to get the results.
    sql_statement = (
        'SELECT * '
        'FROM %s '
        % (mvm.get_table_name(),)
    )
    cursor = connection.cursor()
    cursor.execute(sql_statement)
    raw_result_list = dictfetchall(cursor)

    result_list = [row for row in raw_result_list
            if row['experiment_sample_uid']]

    # We write the following manually specified fields, as well as all the
    # key-value fields.
    csv_field_names = [
        'experiment_sample_uid',
        'experiment_sample_label',
        'uid',
        'position',
        'chromosome',
        'ref',
        'alt',
    ]

    # Add key-value fields
    extra_keys = (
        ref_genome.get_variant_caller_common_map().keys() +
        ref_genome.get_variant_evidence_map().keys()
    )
    csv_field_names.extend(extra_keys)

    # Set structure for checking acceptable fields.
    csv_field_names_set = set(csv_field_names)

    # Write the data.
    writer = csv.DictWriter(csvfile, csv_field_names)
    writer.writeheader()
    for v_idx, variant_data in enumerate(result_list):
        row_data = {
            'experiment_sample_uid': variant_data['experiment_sample_uid'],
            'experiment_sample_label': variant_data['experiment_sample_label'],
            'uid': variant_data['uid'],
            'position': variant_data['position'],
            'chromosome': variant_data['chromosome'],
            'ref': variant_data['ref'],
            'alt': variant_data['alt'],
        }

        # Add key-value data.
        def _add_key_value(data_dict):
            if data_dict is None:
                return
            for key, value in data_dict.iteritems():
                # If this key is not in the header fields, then skip it.
                if not key in csv_field_names_set:
                    continue
                row_data[key] = value
        _add_key_value(variant_data['va_data'])
        _add_key_value(variant_data['vccd_data'])
        _add_key_value(variant_data['ve_data'])

        writer.writerow(row_data)
