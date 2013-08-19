"""
Django sugar for making forms.
"""

from django import forms

from main.models import Project


class ProjectForm(forms.ModelForm):
    """Form used to create a new Project.
    """
    class Meta:
        model = Project
        fields = ('title',)
