from django.shortcuts import render


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


def alignment_list_view(request):
    context = {}
    return render(request, 'alignment_list.html', context)
