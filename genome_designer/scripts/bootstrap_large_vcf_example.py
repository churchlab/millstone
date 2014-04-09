"""
Bootstraps from the fix recoli vcf.
"""

import os
import pickle
import shutil

# Since this script is intended to be used from the terminal, setup the
# environment first so that django and model imports work.
from util import setup_django_env
setup_django_env()

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.model_utils import clean_filesystem_location
from pipeline.variant_effects import get_snpeff_vcf_output_path
from scripts.bootstrap_data import get_or_create_user
from scripts.import_util import import_reference_genome_from_local_file
from settings import PWD as GD_ROOT
from vcf_parser import parse_alignment_group_vcf

FIX_RECOLI_FCF = ''

EXAMPLE_PROJECT_NAME = 'Fix rE. coli'

TEST_DATA_DIR = os.path.join(GD_ROOT, 'test_data')

LARGE_VCF = os.path.join(TEST_DATA_DIR,
        'fix_recoli_variants_snpeff.vcf')

# Uncomment to do subset of ~50 variants.
# LARGE_VCF = os.path.join(TEST_DATA_DIR,
#         'fix_recoli_variants_snpeff_small_test.vcf')

MG1655_REF_GENOME = os.path.join(TEST_DATA_DIR,
        'mg1655_U00096.gbk')

EXPERIMENT_SAMPLE_MODEL_DATA_PICKLE = os.path.join(TEST_DATA_DIR,
        'fix_recoli_variants_snpeff_sample_data.pickle')


def main():
    # Create a User and Project.
    user = get_or_create_user()
    test_project = Project.objects.create(
            title=EXAMPLE_PROJECT_NAME, owner=user.get_profile())
    ref_genome = import_reference_genome_from_local_file(test_project,
            'mg1655', MG1655_REF_GENOME, 'genbank', move=False)

    # Create alignment group and and relate the vcf Dataset to it.
    alignment_group = AlignmentGroup.objects.create(
            label='Fix Recoli Alignment',
            reference_genome=ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)
    vcf_output_path = get_snpeff_vcf_output_path(alignment_group,
            Dataset.TYPE.BWA_ALIGN)
    shutil.copy(LARGE_VCF, vcf_output_path)
    dataset = Dataset.objects.create(
            type=Dataset.TYPE.VCF_FREEBAYES_SNPEFF,
            label=Dataset.TYPE.VCF_FREEBAYES_SNPEFF,
            filesystem_location=clean_filesystem_location(vcf_output_path),
    )
    alignment_group.dataset_set.add(dataset)

    # Import ExperimentSampleo objects, setting specific uid to match
    # the vcf file.
    with open(EXPERIMENT_SAMPLE_MODEL_DATA_PICKLE) as sample_data_fh:
        es_data = pickle.load(sample_data_fh)
        for es in es_data:
            es_obj = ExperimentSample.objects.create(
                uid=es.uid,
                project=test_project,
                label=es.label
            )

            es_obj.data.update(
                {'group':es.group,
                 'well':es.well,
                 'num_reads':es.num_reads})
            es_obj.save()

    parse_alignment_group_vcf(alignment_group,
            Dataset.TYPE.VCF_FREEBAYES_SNPEFF)


if __name__ == '__main__':
    main()
