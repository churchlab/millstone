import json

from django.contrib.auth.decorators import login_required
from django.shortcuts import render

from main.models import AlignmentGroup
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant

from scripts.import_util import import_reference_genome_from_local_file
from scripts.import_util import import_samples_from_targets_file

def home_view(request):
    """The main landing page.
    """
    context = {}
    return render(request, 'home.html', context)


@login_required
def project_list_view(request):
    """The list of projects.
    """
    # NOTE: The 'project_list' template variable is provided via
    # the custom context processor main.context_processors.common_data.
    context = {}
    return render(request, 'project_list.html', context)


@login_required
def project_view(request, project_uid):
    """Overview of a single project.
    """
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project,
    }
    return render(request, 'project.html', context)


@login_required
def reference_genome_list_view(request, project_uid):
    """Shows the ReferenceGenomes and handles creating new
    ReferenceGenomes when requested.
    """
    project = Project.objects.get(uid=project_uid)
    error_string = None

    # If a POST, then we are creating a new genome.
    if request.method == 'POST':
        # TODO: Add more inforative error handling.
        try:
            import_reference_genome_from_local_file(
                    project,
                    request.POST['refGenomeLabel'],
                    request.POST['refGenomeFileLocation'],
                    request.POST['importFileFormat'])
        except Exception as e:
            error_string = 'Import error: ' + str(e)

    # Grab all the ReferenceGenomes for this project.
    ref_genome_list = ReferenceGenome.objects.filter(project=project)


    # Adapt the backend objects to the frontend format.
    fe_ref_genome_list = [{
        'label': obj.label,
        'num_chromosomes': obj.num_chromosomes,
        'total_size': obj.num_bases,
        'annotated': False,
        'parents': [],
        'children': [],
    } for obj in ref_genome_list]

    # HACK: Add some fake data for now.
    fe_ref_genome_list.append({
        'label': 'MG1655-fake',
        'num_chromosomes': '1',
        'total_size': '3.6 Mbp',
        'annotated': False,
        'parents': [],
        'children': ['C321D-fake', 'another'],
    })

    fe_ref_genome_list.append({
        'label': 'C321D-fake',
        'num_chromosomes': '1',
        'total_size': '3.6 Mbp',
        'annotated': True,
        'parents': ['MG1655-fake'],
        'children': [],
    })

    # (Re-)Render the page.
    context = {
        'project': project,
        'ref_genome_list': fe_ref_genome_list,
        'error_string': error_string
    }
    return render(request, 'reference_genome_list.html', context)


@login_required
def sample_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    error_string = None

    # If a POST, then we are creating a new genome.
    if request.method == 'POST':
        # TODO: Add more inforative error handling
        print "Recieved POST!"
        try:
            import_samples_from_targets_file(
                    project,
                    request.FILES['targetsFile'])
        except Exception as e:
            error_string = 'Import error: ' + str(e)
    
    context = {
        'project': project,
        'error_string': error_string
    }
    return render(request, 'sample_list.html', context)



@login_required
def sample_list_targets_template(request):
    """Let the user download a blank sample targets template as a tab
    separated values file (.tsv) so they can fill it in and upload
    it back to the server. 
    """
    context = {}
    return render(request, 'sample_list_targets_template.tsv', context, 
            content_type='text/tab-separated-values')


@login_required
def alignment_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    alignment_list = AlignmentGroup.objects.filter(
            reference_genome__project=project)

    # Adapt the backend objects to the frontend format.
    fe_alignment_list = [{
        'label': alignment_group_obj.label,
        'reference_genome': alignment_group_obj.reference_genome,
        'sample_desc': str(len(alignment_group_obj.experimentsampletoalignment_set.all())) + ' samples',
        'start_time': alignment_group_obj.start_time,
        'end_time': alignment_group_obj.end_time,
    } for alignment_group_obj in alignment_list]

    context = {
        'project': project,
        'alignment_list': fe_alignment_list
    }
    return render(request, 'alignment_list.html', context)


@login_required
def variant_set_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project
    }
    return render(request, 'variant_set_list.html', context)


@login_required
def variant_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)

    # Fetch the list of variants and render it into the dom as json.
    # The data will be displayed to the user via the javascript DataTables
    # component.
    variant_list = Variant.objects.all()

    fe_variant_list = [{
        'type': obj.type,
        'chromosome': obj.chromosome,
        'position': obj.position,
        'ref_value': obj.ref_value,
        'alt_value': obj.alt_value,
    } for obj in variant_list]

    # Configuration for display using jquery.datatables.js.
    variant_field_config = [
        {'mData': 'type', 'sTitle': 'Type'},
        {'mData': 'chromosome', 'sTitle': 'Chromosome'},
        {'mData': 'position', 'sTitle': 'Position'},
        {'mData': 'ref_value', 'sTitle': 'Ref Value'},
        {'mData': 'alt_value', 'sTitle': 'Alt Value'}
    ]

    variant_list_json = json.dumps({
        'variant_list': fe_variant_list,
        'field_config': variant_field_config
    })

    context = {
        'project': project,
        'variant_list_json': variant_list_json
    }
    return render(request, 'variant_list.html', context)


@login_required
def gene_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project
    }
    return render(request, 'gene_list.html', context)


@login_required
def goterm_list_view(request, project_uid):
    project = Project.objects.get(uid=project_uid)
    context = {
        'project': project
    }
    return render(request, 'goterm_list.html', context)
