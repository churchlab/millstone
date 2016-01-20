from collections import defaultdict
from collections import OrderedDict

from main.models import AlignmentGroup
from main.models import ExperimentSampleToAlignment
from variants.materialized_variant_filter import lookup_variants


def get_tenaillon_de_novo_variants_by_sample():
    alignment_label = 'All Lines to REL1206'
    alignment = AlignmentGroup.objects.get(label=alignment_label)

    return get_de_novo_variants_by_sample(alignment)


def get_de_novo_variants_by_sample(alignment):

    query_args = {
            'filter_string': 'INFO_METHOD=DE_NOVO_ASSEMBLY',
            'melted': False
    }
    de_novo_assembly_variants = lookup_variants(
            query_args, alignment.reference_genome, alignment).result_list

    completed_samples = [sample_alignment.experiment_sample
            for sample_alignment in ExperimentSampleToAlignment.objects.filter(
                    alignment_group=alignment) if (
                    sample_alignment.data.get('assembly_status', None) == (
                    ExperimentSampleToAlignment.ASSEMBLY_STATUS.COMPLETED))]

    variants_by_sample_label = defaultdict(list)
    for variant in de_novo_assembly_variants:
        variants_by_sample_label[variant['EXPERIMENT_SAMPLE_LABEL']].append(
                variant)

    variants_by_completed_sample_label = {}
    for sample in completed_samples:
        variants_by_completed_sample_label[sample.label] = (
                variants_by_sample_label[sample.label])

    return variants_by_completed_sample_label


def readable_variant(variant):
    alt = variant['ALT']
    alt_repr = alt if alt.startswith('<') else (str(len(alt)) + ' Bases')

    ref = variant['REF']
    ref_repr = str(len(ref)) + ' Bases'

    position = variant['POSITION']
    position_repr = str(position)

    return '\t'.join([
            'POSITION:' + position_repr,
            'REF: ' + ref_repr,
            'ALT: ' + alt_repr
    ])


def make_report(sv_by_sample_dict):

    for variant_list in sv_by_sample_dict.values():
        variant_list.sort(key=lambda v: v['POSITION'])

    readable_sv_by_sample_dict = {s: map(readable_variant, v) for s, v in
            sv_by_sample_dict.items()}

    ordered_dict = OrderedDict(sorted(readable_sv_by_sample_dict.items(),
            key=lambda x: x[0]))

    return '\n\n'.join([
        label + '\n\t' + '\n\t'.join(v_list)
        for label, v_list in ordered_dict.items()])
