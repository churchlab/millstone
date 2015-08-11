"""
Template context processors specific to Genome Designer.

These are variables that are automatically made available to the template
renderer.
"""

import re
import subprocess

from main.models import Project
from django.conf import settings


def common_data(request):
    """Returns context variables with data common to every view, for example
    the projects belonging to a user.
    """
    # Data that is accessible from every page.
    data = {
        'FLAG__GENOME_FINISHING_ENABLED':
                settings.FLAG__GENOME_FINISHING_ENABLED,
        'free_disk_space_str': _get_free_disk_space_str()
    }

    # Add data that is user-specific.
    user_specific_data = {}
    if hasattr(request, 'user') and request.user.is_authenticated():
        user_specific_data['project_list'] = Project.objects.filter(
                owner=request.user)

    # Put it all together and ship it.
    data.update(user_specific_data)
    return data


def _get_free_disk_space_str():
    """Returns html string representing free disk space.

    Adds danger highlight if less than 10% remaining.
    """
    # Run native df.
    df = subprocess.Popen(
            ['df', '-h', settings.MEDIA_ROOT], stdout=subprocess.PIPE)
    df_output = df.communicate()[0]

    # Parse output.
    parts = df_output.split('\n')[1].split()
    size = parts[1]
    available = parts[3]

    size_int = int(re.match(r'[0-9]+', size).group())
    available_int = int(re.match(r'[0-9]+', available).group())
    capacity_int = size_int * 100 / available_int

    # Maybe add red error text style to capacity.
    available_str = available
    if capacity_int < 10:
        available_str = (
                '<span class="gd-error-text">' + available_str + '</span>')

    return '{available} out of {size}'.format(
            available=available_str, size=size)


def aws_settings(request):
    return {
        'AWS_SERVER_PUBLIC_KEY': settings.AWS_SERVER_PUBLIC_KEY,
        'S3_BUCKET': settings.S3_BUCKET,
    }


def demo_settings(request):
    return {
        'is_demo': settings.DEMO_MODE
    }
