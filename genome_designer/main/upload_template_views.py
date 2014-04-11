"""
Handlers for fetching upload templates.
"""

from django.shortcuts import render


SAMPLE_LIST_BROWSER_UPLOAD_TEMPLATE = 'sample_list_browser_upload_template.tsv'


def sample_list_targets_template(request):
    """Let the user download a blank sample targets template as a tab
    separated values file (.tsv) so they can fill it in and upload
    it back to the server.
    """
    context = {}
    return render(request, 'sample_list_targets_template.tsv', context,
            content_type='text/tab-separated-values')


def variant_set_upload_template(request):
    """Let the user download a blank variant set template as a blank
    VCF file to be filled in.
    """
    context = {}
    return render(request, 'variant_set_upload_template.vcf', context,
            content_type='text/tab-separated-values')


def sample_list_browser_upload_template(request):
    """Template that allows the user to indicate the names of samples they
    will updload through the browser form.
    """
    context = {}
    return render(request, SAMPLE_LIST_BROWSER_UPLOAD_TEMPLATE, context,
            content_type='text/tab-separated-values')
