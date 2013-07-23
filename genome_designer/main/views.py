from django.shortcuts import render

from main.models import AlignmentGroup


def home_view(request):
    """The main landing page.
    """
    context = {}
    return render(request, 'home.html', context)


def project_list_view(request):
    """The list of projects.
    """
    context = {}
    return render(request, 'project_list.html', context)


def reference_genome_list_view(request):
    context = {}
    return render(request, 'reference_genome_list.html', context)


def sample_list_view(request):
    context = {}
    return render(request, 'sample_list.html', context)
    
def sample_list_upload_template(request):
    """Let the user download a blank sample upload template as a tab
    separated values file (.tsv) and allow them to fill it in and upload
    it back to the server. 
    """
    
    context = {}
    return render(request, 'sample_list_upload_template.tsv', context, 
            content_type='text/tab-separated-values')


def alignment_list_view(request):
    alignment_list = AlignmentGroup.objects.all()

    # Adapt the backend objects to the frontend format.
    fe_alignment_list = [{
        'label': alignment_group_obj.label,
        'reference_genome': alignment_group_obj.reference_genome,
        'sample_desc': str(len(alignment_group_obj.experimentsampletoalignment_set.all())) + ' samples',
        'start_time': alignment_group_obj.start_time,
        'end_time': alignment_group_obj.end_time,
    } for alignment_group_obj in alignment_list]

    context = {
        'alignment_list': fe_alignment_list
    }
    return render(request, 'alignment_list.html', context)


def variant_set_list_view(request):
    context = {}
    return render(request, 'variant_set_list.html', context)


def variant_list_view(request):
    context = {}
    return render(request, 'variant_list.html', context)


def gene_list_view(request):
    context = {}
    return render(request, 'gene_list.html', context)


def goterm_list_view(request):
    context = {}
    return render(request, 'goterm_list.html', context)
