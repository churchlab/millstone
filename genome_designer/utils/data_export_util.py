"""
Methods for exporting data.
"""

import csv
from datetime import datetime
import os
import re
import StringIO
import zipfile

from Bio import SeqIO
from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from utils import lowercase_underscore
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
EXPORT_CONTIGS_TEMPLATE_PATH = os.path.join(settings.PWD, 'test_data',
        'export_contigs_vcf_template.vcf')

# We use a fake single sample. The reason for this is that the external library
# reference_genome_maker expects there to be a sample in case there are
# multiple alts. We're not necessarily set on continuing to do it this way,
# but we're happy enough with it for now.
PLACEHOLDER_SAMPLE_NAME = 'fake_sid'


def export_contig_list_as_vcf(contig_list, vcf_dest_path_or_filehandle):
    """Exports a list of contigs as a vcf
    """
    # Allow dest input as path or filehandle.
    if (isinstance(vcf_dest_path_or_filehandle, str) or
            isinstance(vcf_dest_path_or_filehandle, unicode)):
        out_vcf_fh = open(vcf_dest_path_or_filehandle, 'w')
    else:
        out_vcf_fh = vcf_dest_path_or_filehandle

    # Assert all contigs come from same sample
    assert len(set([c.experiment_sample_to_alignment.experiment_sample.uid
            for c in contig_list])) == 1

    # The vcf.Writer() requires a template vcf in order to structure itself.
    # We read in a generic template to satisfy this.
    with open(EXPORT_CONTIGS_TEMPLATE_PATH) as template_fh:
        vcf_template = vcf.Reader(template_fh)

    # Set samples and sample_indexes for template
    contig_0 = contig_list[0]
    sample_uid = (
            contig_0.experiment_sample_to_alignment.experiment_sample.uid)

    vcf_template.samples = [sample_uid]
    sample_indexes = {sample_uid: 0}


    vcf_writer = vcf.Writer(out_vcf_fh, vcf_template)
    for contig in contig_list:
        assert contig.chromosome
        assert contig.reference_insertion_endpoints
        assert contig.contig_insertion_endpoints

        ref_left, ref_right = contig.reference_insertion_endpoints
        contig_left, contig_right = contig.contig_insertion_endpoints

        # Get Seqrecord
        contig_fasta = get_dataset_with_type(
                contig,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
        with open(contig_fasta) as fh:
                contig_seqrecord = SeqIO.parse(fh, 'fasta').next()

        # Determine whether contig is reverse complement relative to reference
        is_reverse = contig.metadata.get('is_reverse', False)

        # Extract cassette sequence from contig
        if is_reverse:
            cassette_sequence = str(contig_seqrecord.seq.reverse_complement()[
                contig_left:contig_right])
        else:
            cassette_sequence = str(contig_seqrecord.seq[
                    contig_left:contig_right])

        if ref_left > ref_right:
            bases_to_peel_back = ref_left - ref_right
            ref_genome_fasta = get_dataset_with_type(
                    contig.parent_reference_genome,
                    Dataset.TYPE.REFERENCE_GENOME_FASTA
                            ).get_absolute_location()
            with open(ref_genome_fasta) as fh:
                ref_seqrecord_iter = SeqIO.parse(fh, 'fasta')
                ref_seqrecord = None
                for seqrecord in ref_seqrecord_iter:
                    if seqrecord.id == contig.chromosome:
                        ref_seqrecord = seqrecord
                        break
                assert ref_seqrecord is not None

            peel_back_sequence = str(ref_seqrecord.seq[
                    (ref_left - bases_to_peel_back):ref_left])

            pos = ref_left - bases_to_peel_back + 1
            ref_value = ''
            alt_value = peel_back_sequence + cassette_sequence

        elif ref_right > ref_left:
            ref_genome_fasta = get_dataset_with_type(
                    contig.parent_reference_genome,
                    Dataset.TYPE.REFERENCE_GENOME_FASTA
                            ).get_absolute_location()
            with open(ref_genome_fasta) as fh:
                ref_seqrecord_iter = SeqIO.parse(fh, 'fasta')
                ref_seqrecord = None
                for seqrecord in ref_seqrecord_iter:
                    if seqrecord.id == contig.chromosome:
                        ref_seqrecord = seqrecord
                        break
                assert ref_seqrecord is not None

            # + 1 to insert AFTER end of reference
            pos = ref_left + 1
            ref_value = str(ref_seqrecord.seq[
                    ref_left + 1: ref_right + 1])
            alt_value = cassette_sequence

        elif ref_right == ref_left:
            pos = ref_left + 1
            ref_value = ''
            alt_value = cassette_sequence

        # In case of no alt this contig represents a deletion which
        # is represented by a '<DEL>' alt field in vcf format
        alt_value = alt_value if alt_value else '<DEL>'

        record = vcf.model._Record(
                contig.chromosome,
                pos,
                contig.uid,
                ref_value,
                (alt_value,),
                1, # QUAL
                [], # FILTER
                {'contig_uid': contig.uid}, # INFO
                'GT:DP:RO:QR:AO:QA:GL', # FORMAT
                sample_indexes, # sample_indexes
        )
        # Add a placeholder sample.
        calldata_type = vcf.model.make_calldata_tuple(['GT'])
        placeholder_sample_data = calldata_type(GT='1/1')
        record.samples = [vcf.model._Call(
                record, sample_uid, placeholder_sample_data)]
        vcf_writer.write_record(record)


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


def export_project_as_zip(project):
    """Compresses project and generates link for download.
    """
    # Common root for exported zip file.
    COMMON_EXPORT_ROOT = 'millstone_export'

    # Delete previous export files to avoid overwhelming space.
    for f in os.listdir(settings.TEMP_FILE_ROOT):
        if re.match(COMMON_EXPORT_ROOT, f):
            os.remove(os.path.join(settings.TEMP_FILE_ROOT, f))

    project_export_zip_name = (
            '{common_root}_{proj_title}_{timestamp}.zip'.format(
                    common_root=COMMON_EXPORT_ROOT,
                    proj_title=lowercase_underscore(project.title[:20]),
                    timestamp=datetime.now().strftime('%Y_%m_%d_%H%M')))
    project_zip_dest = os.path.join(
            settings.TEMP_FILE_ROOT, project_export_zip_name)
    with zipfile.ZipFile(project_zip_dest, 'w', allowZip64=True) as ziph:
        project_root_dir = project.get_model_data_dir()
        for root, dirs, files in os.walk(project_root_dir):
            for file in files:
                full_path = os.path.join(root, file)
                # Write to destination in archive that omits part of path
                # before project uid.
                archive_target_path = (
                        re.search(project.uid + '.*', full_path).group())
                ziph.write(full_path, archive_target_path)

    return os.path.join('/tmp', project_export_zip_name)
