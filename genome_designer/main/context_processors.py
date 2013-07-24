"""
Template context processors specific to Genome Designer.

These are variables that are automatically made available to the template
renderer.
"""

from main.models import Project

def common_data(request):
    """Returns context variables with data common to every view, for example
    the projects belonging to a user.
    """
    if not hasattr(request, 'user') or not request.user.is_authenticated():
        return {}
    return {
        'project_list': Project.objects.filter(owner=request.user)
    }
