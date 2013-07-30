from django.conf.urls.defaults import include
from django.conf.urls.defaults import patterns
from django.conf.urls.defaults import url

urlpatterns = patterns('',
    url(r'^$', 'genome_designer.main.views.home_view'),

    # Project-specific views
    url(r'^projects$',
            'genome_designer.main.views.project_list_view'),
    url(r'^projects/([\w-]+)$',
            'genome_designer.main.views.project_view'),

    # Reference genomes
    url(r'^projects/([\w-]+)/refgenomes$',
            'genome_designer.main.views.reference_genome_list_view'),
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)$',
            'genome_designer.main.views.reference_genome_view'),

    # Alignments
    url(r'^projects/([\w-]+)/alignments$',
            'genome_designer.main.views.alignment_list_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)$',
            'genome_designer.main.views.alignment_view'),
    url(r'^projects/([\w-]+)/alignments/create$',
            'genome_designer.main.views.alignment_create_view'),

    # Variant sets
    url(r'^projects/([\w-]+)/sets$',
            'genome_designer.main.views.variant_set_list_view'),

    # Samples
    url(r'^projects/([\w-]+)/samples$',
            'genome_designer.main.views.sample_list_view'),

    # Variants
    url(r'^projects/([\w-]+)/variants$',
            'genome_designer.main.views.variant_list_view'),

    # Genes
    url(r'^projects/([\w-]+)/genes$',
            'genome_designer.main.views.gene_list_view'),

    # GO terms
    url(r'^projects/([\w-]+)/goterms$',
            'genome_designer.main.views.goterm_list_view'),


    ############################################################################
    # Templates
    ############################################################################

    url(r'^templates/sample_list_targets_template.tsv$',
            'genome_designer.main.views.sample_list_targets_template'),


    ############################################################################
    # Auth
    ############################################################################

    # django-registration defaults (further delgates to django.contrib.auth.url)
    (r'^accounts/', include('registration.backends.simple.urls')),
)
