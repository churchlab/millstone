"""
Methods for exporting data.
"""

import csv
import os
import StringIO

from django.conf import settings
import vcf

from variants.materialized_variant_filter import get_variants_that_pass_filter
from variants.materialized_view_manager import MeltedVariantMaterializedViewManager


CORE_VARIANT_KEYS = [
    'UID',
    'POSITION',
    'CHROMOSOME',
    'REF',
    'ALT',
    'EXPERIMENT_SAMPLE_LABEL',
    'VARIANT_SET_LABEL'
]


def export_melted_variant_view(alignment_group, filter_string):
    """Generator that yields rows of a csv file.

    Args:
        ref_genome: ReferenceGenome these Variants belong to.
        filter_string: Limit the returned Variants to those that match this
            filter.
    """
    ref_genome = alignment_group.reference_genome

    mvm = MeltedVariantMaterializedViewManager(ref_genome)
    mvm.create_if_not_exists_or_invalid()

    # We'll perform a query, using any filter_string provided.
    query_args = {}
    query_args['filter_string'] = filter_string
    query_args['select_all'] = True
    query_args['act_as_generator'] = True
    variant_iterator = get_variants_that_pass_filter(
            query_args, ref_genome, alignment_group=alignment_group)

    # We write the core keys and key-values specific to this ref_genome.
    ref_genome_specific_data_keys = (
        ref_genome.get_variant_caller_common_map().keys() +
        ref_genome.get_variant_evidence_map().keys() +
        ref_genome.get_variant_alternate_map().keys() +
        ref_genome.get_experiment_sample_map().keys()
    )

    all_keys = CORE_VARIANT_KEYS + ref_genome_specific_data_keys

    # Our convention is to capitalize all keys.
    csv_field_names = [key.upper() for key in all_keys]

    # Create a set that is used to quickly check which fields to return.
    csv_field_names_set = set(csv_field_names)

    # Create a csv writer that uses a StringIO buffer.
    # Every time we write a row, we flush the buffer and yield the data.
    output_buffer = StringIO.StringIO()
    writer = csv.DictWriter(output_buffer, csv_field_names)

    # Write header
    writer.writeheader()
    output_buffer.seek(0)
    data = output_buffer.read()
    output_buffer.truncate(0)
    yield data

    for variant_data in variant_iterator:
        row_data = {}
        for key in CORE_VARIANT_KEYS:
            row_data[key] = variant_data[key]

        # Add key-value data.
        def _add_key_value(data_dict):
            if data_dict is None:
                return
            for key, value in data_dict.iteritems():
                # If this key is not in the header fields, then skip it.
                if not key in csv_field_names_set:
                    continue
                row_data[key] = value
        _add_key_value(variant_data['VA_DATA'])
        _add_key_value(variant_data['VCCD_DATA'])
        _add_key_value(variant_data['VE_DATA'])

        writer.writerow(row_data)
        output_buffer.seek(0)
        data = output_buffer.read()
        output_buffer.truncate(0)
        yield data


VCF_TEMPLATE_PATH = os.path.join(settings.PWD, 'test_data', 'vcf_template.vcf')

# We use a fake single sample. The reason for this is that the external library
# reference_genome_maker expects there to be a sample in case there are
# multiple alts. We're not necessarily set on continuing to do it this way,
# but we're happy enough with it for now.
PLACEHOLDER_SAMPLE_NAME = 'fake_sid'


def export_variant_set_as_vcf(variant_set, vcf_dest_path_or_filehandle):
    """Exports a VariantSet as a vcf.
    """
    # Allow dest input as path or filehandle.
    if isinstance(vcf_dest_path_or_filehandle, str):
        out_vcf_fh = open(vcf_dest_path_or_filehandle, 'w')
    else:
        out_vcf_fh = vcf_dest_path_or_filehandle

    # The vcf.Writer() requires a template vcf in order to structure itself.
    # We read in a generic template to satisfy this.
    with open(VCF_TEMPLATE_PATH) as template_fh:
        vcf_template = vcf.Reader(template_fh)

    vcf_writer = vcf.Writer(out_vcf_fh, vcf_template)
    for variant in variant_set.variants.all():
        alts = variant.variantalternate_set.all()
        assert len(alts) == 1, "Only support variants with exactly one alt."
        alt_value = alts[0].alt_value

        record = vcf.model._Record(
                variant.chromosome.seqrecord_id,
                variant.position,
                variant.uid,
                variant.ref_value,
                (alt_value,),
                1, # QUAL
                [], # FILTER
                {}, # INFO
                'GT:DP:RO:QR:AO:QA:GL', # FORMAT
                {PLACEHOLDER_SAMPLE_NAME: 0}, # sample_indexes
        )
        # Add a placeholder sample.
        calldata_type = vcf.model.make_calldata_tuple(['GT'])
        placeholder_sample_data = calldata_type(GT='1/1')
        record.samples = [vcf.model._Call(
                record, PLACEHOLDER_SAMPLE_NAME, placeholder_sample_data)]
        vcf_writer.write_record(record)
