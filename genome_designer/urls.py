from django.conf.urls import include
from django.conf.urls import patterns
from django.conf.urls import url
from django.views.generic import RedirectView

import settings

urlpatterns = patterns('',
    url(r'^$', 'genome_designer.main.views.home_view'),

    # Project-specific views
    url(r'^projects$',
            'genome_designer.main.views.project_list_view'),
    url(r'^projects/create$',
            'genome_designer.main.views.project_create_view'),
    url(r'^projects/([\w-]+)$',
            'genome_designer.main.views.project_view'),
    url(r'^projects/([\w-]+)/delete$',
            'genome_designer.main.views.project_delete'),


    # Tab base views.
    url(r'^projects/([\w-]+)/data$',
            'genome_designer.main.views.project_view'),
    url(r'^projects/([\w-]+)/align$',
            'genome_designer.main.views.tab_root_align'),
    url(r'^projects/([\w-]+)/analyze$',
            'genome_designer.main.views.tab_root_analyze'),


    # Reference genomes
    url(r'^projects/([\w-]+)/refgenomes$',
            'genome_designer.main.views.reference_genome_list_view'),
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)$',
            'genome_designer.main.views.reference_genome_view'),

    # Alignments
    url(r'^projects/([\w-]+)/alignments$',
            'genome_designer.main.views.alignment_list_view'),
    url(r'^projects/([\w-]+)/alignments/create$',
            'genome_designer.main.views.alignment_create_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)$',
            'genome_designer.main.views.alignment_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)/samplealign/([\w-]+)/error$',
            'genome_designer.main.views.sample_alignment_error_view'),


    # Variant sets
    url(r'^projects/([\w-]+)/sets$',
            'genome_designer.main.views.variant_set_list_view'),
    url(r'^projects/([\w-]+)/sets/([\w-]+)$',
            'genome_designer.main.views.variant_set_view'),


    # Samples
    url(r'^projects/([\w-]+)/samples$',
            'genome_designer.main.views.sample_list_view'),

    # Variants
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)/variants/([\w-]+)$',
            'genome_designer.main.views.single_variant_view'),


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

    url(r'^templates/variant_set_upload_template.vcf$',
            'genome_designer.main.views.variant_set_upload_template'),


    ############################################################################
    # Auth
    ############################################################################

    # django-registration defaults (further delgates to django.contrib.auth.url)
    (r'^accounts/', include('registration.backends.simple.urls')),

    # The default behavior of registration is redirect to 'users/<username>'.
    # For now let's catch this request here and just redirect to '/'.
    (r'^users/', RedirectView.as_view(url='/')),


    ############################################################################
    # XHR Actions
    ############################################################################

    url(r'^_/sets/exportcsv$',
            'genome_designer.main.xhr_handlers.export_variant_set_as_csv'),
    url(r'^_/variants$',
            'genome_designer.main.xhr_handlers.get_variant_list'),
    url(r'^_/variants/modify_set_membership$',
            'genome_designer.main.xhr_handlers.modify_variant_in_set_membership'),
)

if settings.DEBUG:
    from django.conf.urls.static import static
    urlpatterns += static('jbrowse', document_root=settings.JBROWSE_ROOT)
