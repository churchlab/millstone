from django.conf.urls.defaults import patterns
from django.conf.urls.defaults import url

urlpatterns = patterns('',
    url(r'^$', 'genome_designer.main.views.home_view'),
    url(r'^projects$', 'genome_designer.main.views.project_list_view'),
    url(r'^refgenomes$', 'genome_designer.main.views.reference_genome_list_view'),
    url(r'^alignments$', 'genome_designer.main.views.alignment_list_view'),
    url(r'^sets$', 'genome_designer.main.views.variant_set_list_view'),
    url(r'^samples$', 'genome_designer.main.views.sample_list_view'),
    url(r'^samples/sample_list_upload_template.tsv$', 'genome_designer.main.views.sample_list_upload_template'),
    url(r'^variants$', 'genome_designer.main.views.variant_list_view'),
    url(r'^genes$', 'genome_designer.main.views.gene_list_view'),
    url(r'^goterms$', 'genome_designer.main.views.goterm_list_view'),
)
