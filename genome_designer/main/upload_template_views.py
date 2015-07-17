"""
Handlers for fetching upload templates.
"""

from django.shortcuts import render


CONTENT_TYPE__CSV = 'text/csv'

TEMPLATE__SAMPLES_BATCH_IMPORT_FROM_SERVER = 'sample_list_targets_template.csv'

TEMPLATE__SAMPLES_BATCH_BROWSER_UPLOAD = (
        'sample_list_browser_upload_template.csv')


def sample_list_targets_template(request):
    """Let the user download a blank sample targets template as a
    comma-separated values file (.csv) so they can fill it in and upload
    it back to the server.
    """
    context = {}
    return render(request, TEMPLATE__SAMPLES_BATCH_IMPORT_FROM_SERVER, context,
            content_type=CONTENT_TYPE__CSV)


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
    return render(request, TEMPLATE__SAMPLES_BATCH_BROWSER_UPLOAD, context,
            content_type=CONTENT_TYPE__CSV)
