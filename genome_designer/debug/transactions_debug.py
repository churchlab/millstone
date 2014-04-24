"""
Experimenting with Django transactions.
"""

import os
import random
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from django.contrib.auth.models import User
from django.db import transaction

from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from scripts.bootstrap_data import REF_GENOME_1_LABEL
from scripts.bootstrap_data import TEST_PROJECT_NAME
from scripts.bootstrap_data import TEST_USERNAME

def main():
    try:
        user = User.objects.get(username=TEST_USERNAME)
    except User.DoesNotExist:
        user = User.objects.create_user(
                TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

    ref_genome, ref_genome_created = ReferenceGenome.objects.get_or_create(
            label=REF_GENOME_1_LABEL)

    test_project = Project.objects.create(
            title='deleteme', owner=user.get_profile())

    # var_list = []
    # for pos in range(3):
    #     var_list.append(Variant(
    #             type=Variant.TYPE.TRANSITION,
    #             reference_genome=ref_genome,
    #             chromosome='chrom',
    #             position=pos,
    #             ref_value='A'))
    # Variant.objects.bulk_create(var_list)

    # with transaction.commit_on_success():
    #     for pos in range(3):
    #         Variant.objects.create(
    #                 type=Variant.TYPE.TRANSITION,
    #                 reference_genome=ref_genome,
    #                 chromosome='chrom',
    #                 position=pos,
    #                 ref_value='A')

if __name__ == '__main__':
    main()
