"""
Convenience module for debugging scripts.
"""

# Since this script is intended to be used from the terminal, setup the
# environment first so that django and model imports work.
from util import setup_django_env
setup_django_env()

from main.models import AlignmentGroup
from main.models import ExperimentSampleToAlignment
from scripts.snp_callers import run_snp_calling_pipeline_internal
from scripts.alignment_pipeline import align_with_bwa

def debug_snp_callers():
    alignment_group = AlignmentGroup.objects.get(
            uid='54305db0')

    run_snp_calling_pipeline_internal(alignment_group)


def debug_align_with_bwa():
    alignment_group = AlignmentGroup.objects.get(
            uid='54305db0')

    sample_alignment = ExperimentSampleToAlignment.objects.get(uid='7ab59bc7')

    align_with_bwa(alignment_group, sample_alignment=sample_alignment)


if __name__ == '__main__':
    debug_align_with_bwa()
