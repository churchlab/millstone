from Bio import SeqIO
import pyinter
import pickle

from main.model_utils import get_dataset_with_type

GBK_FEATURES_TO_EXTRACT = [
        'CDS',
        'mobile_element',
        'tRNA',
        'rRNA']

def get_genbank_features_with_type(genbank_path, feature_type):

    chrom_intervals = {}
    with open(genbank_path, 'r') as fh:
        for seq_record in SeqIO.parse(fh, 'genbank'):
            interval_list = []
            for f in seq_record.features:
                if f.type == feature_type:
                    interval_list.append((f, f.extract(seq_record.seq)))

            chrom_intervals[seq_record.id] = interval_list

    return chrom_intervals


def generate_genbank_mobile_element_multifasta(genbank_path, output_fasta_path):
    """Extract mobile elements from a genbank-annotated genome into a multifasta
    record, for use with SV and mobile element-calling code.
    """
    me_features = get_genbank_features_with_type(
            genbank_path, 'mobile_element')
    me_sequences = {}
    for chrom, feature_seq_tuples in me_features.items():
        for feature, seq in feature_seq_tuples:
            if 'mobile_element_type' not in feature.qualifiers:
                continue
            me_type = feature.qualifiers['mobile_element_type'][0]
            me_type = me_type.replace(' ', '_')
            if me_type not in me_sequences:
                me_sequences[me_type] = seq

    with open(output_fasta_path, 'w') as fh:
        for me_type, seq in me_sequences.items():
            fh.write('>ME_%s\n' % me_type)
            fh.write(str(seq))
            fh.write('\n')


def generate_gbk_feature_index(genbank_path, feature_index_output_path):
    """
    Create a pickled pyinterval index of genbank features so we can pull
    them quickly.
    """

    gbk_feature_list = []
    with open(genbank_path, 'r') as fh:
        for seq_record in SeqIO.parse(fh, 'genbank'):
            interval_list = []

            for f in seq_record.features:

                if f.type not in GBK_FEATURES_TO_EXTRACT:
                    continue

                f_ivl = pyinter.closedopen(
                        f.location.start, f.location.end)
                f_ivl.type = f.type

                if 'gene' in f.qualifiers:
                    f_ivl.name = f.qualifiers['gene'][0]
                elif 'mobile_element_type' in f.qualifiers:
                    f_ivl.name = f.qualifiers['mobile_element_type'][0]
                # For now, if the gene has no '.name' or '.mobile_element_type', ignore
                gbk_feature_list.append(f_ivl)

    with open(feature_index_output_path, 'w') as fh:
        pickle.dump(gbk_feature_list, fh)

