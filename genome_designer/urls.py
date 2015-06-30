from django.conf import settings
from django.conf.urls import include
from django.conf.urls import patterns
from django.conf.urls import url
from django.views.generic import RedirectView
from django.views.generic.base import TemplateView

from main.views import RegistrationViewWrapper

urlpatterns = patterns('',
    url(r'^$', 'main.views.home_view'),

    # Project-specific views
    url(r'^projects$',
            'main.views.project_list_view'),
    url(r'^projects/create$',
            'main.views.project_create_view'),
    url(r'^projects/([\w-]+)$',
            'main.views.project_view'),

    # Tab base views.
    url(r'^projects/([\w-]+)/data$',
            'main.views.project_view'),
    url(r'^projects/([\w-]+)/analyze$',
            'main.views.tab_root_analyze'),
    url(r'^projects/([\w-]+)/analyze/([\w-]+)$',
            'main.views.tab_root_analyze'),
    url(r'^projects/([\w-]+)/analyze/([\w-]+)/([\w-]+)$',
            'main.views.tab_root_analyze'),

    # Genome finishing
    url(r'^projects/([\w-]+)/genomefinishing$',
            'main.views.genome_finish_view'),

    # Reference genomes
    url(r'^projects/([\w-]+)/refgenomes$',
            'main.views.reference_genome_list_view'),
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)$',
            'main.views.reference_genome_view'),

    # Alignments
    url(r'^projects/([\w-]+)/alignments$',
            'main.views.alignment_list_view'),
    url(r'^projects/([\w-]+)/alignments/create$',
            'main.views.alignment_create_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)$',
            'main.views.alignment_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)/samplealign/([\w-]+)/error$',
            'main.views.sample_alignment_error_view'),

    # Variant sets
    url(r'^projects/([\w-]+)/sets$',
            'main.views.variant_set_list_view'),
    url(r'^projects/([\w-]+)/sets/([\w-]+)$',
            'main.views.variant_set_view'),

    # Samples
    url(r'^projects/([\w-]+)/samples$',
            'main.views.sample_list_view'),
    url(r'^projects/([\w-]+)/samples/([\w-]+)/fastqc/read([\d+])$',
            'main.views.fastqc_view'),


    # Demo Splash
    url(r'^demo_splash$',
            'main.views.demo_splash_view'),


    ###########################################################################
    # Templates
    ###########################################################################

    url(r'^templates/sample_list_browser_upload_template.tsv$',
            'main.upload_template_views.sample_list_browser_upload_template'),

    url(r'^templates/sample_list_targets_template.tsv$',
            'main.upload_template_views.sample_list_targets_template'),

    url(r'^templates/variant_set_upload_template.vcf$',
            'main.upload_template_views.variant_set_upload_template'),


    ###########################################################################
    # XHR Actions
    ###########################################################################

    url(r'^_/sets$',
            'main.xhr_handlers.get_variant_set_list'),
    url(r'^_/sets/create$',
            'main.xhr_handlers.create_variant_set'),
    url(r'^_/sets/print_mage_oligos$',
            'main.xhr_handlers.print_mage_oligos_for_variant_set'),
    url(r'^_/sets/generate_new_ref_genome$',
            'main.xhr_handlers.generate_new_ref_genome_for_variant_set'),

    url(r'^_/variants$',
            'main.xhr_handlers.get_variant_list'),
    url(r'^_/variants/modify_set_membership$',
            'main.xhr_handlers.modify_variant_in_set_membership'),
    url(r'^_/variants/refresh_materialized_variant_table$',
            'main.xhr_handlers.refresh_materialized_variant_table'),
    url(r'^_/variants/export_as_csv$',
            'main.xhr_handlers.export_variants_as_csv'),
    url(r'^_/variants/is_materialized_view_valid$',
            'main.xhr_handlers.is_materialized_view_valid'),
    url(r'^_/variants/delete$',
        'main.xhr_handlers.variant_sets_delete'),
    url(r'^_/variants/save_filter$',
        'main.xhr_handlers.save_variant_filter'),
    url(r'^_/variants/delete_filter$',
        'main.xhr_handlers.delete_variant_filter'),

    url(r'^_/alignmentgroups$',
            'main.xhr_handlers.get_alignment_groups'),
    url(r'^_/alignmentgroups/delete$',
            'main.xhr_handlers.alignment_groups_delete'),

    url(r'^_/alignments/generate_contigs$',
            'main.xhr_handlers.generate_contigs'),

    url(r'^_/samples$',
            'main.xhr_handlers.get_samples'),
    url(r'^_/samples/upload_single_sample$',
            'main.xhr_handlers.upload_single_sample'),
    url(r'^_/samples/create_from_server_location$',
            'main.xhr_handlers.create_samples_from_server_location'),
    url(r'^_/samples/samples_upload_through_browser_template$',
            'main.xhr_handlers.samples_upload_through_browser_template'),
    url(r'^_/samples/samples_upload_through_browser_sample_data$',
            'main.xhr_handlers.samples_upload_through_browser_sample_data'),
    url(r'^_/samples/get_samples_awaiting_upload$',
            'main.xhr_handlers.get_samples_awaiting_upload'),
    url(r'^_/samples/delete$',
            'main.xhr_handlers.samples_delete'),

    url(r'^_/genes$',
            'main.xhr_handlers.get_gene_list'),

    url(r'^projects/([\w-]+)/delete$',
            'main.xhr_handlers.project_delete'),

    url(r'^_/ref_genomes$',
            'main.xhr_handlers.get_ref_genomes'),
    url(r'^_/ref_genomes/upload_through_browser$',
            'main.xhr_handlers.create_ref_genome_from_browser_upload'),
    url(r'^_/ref_genomes/create_from_server_location$',
            'main.xhr_handlers.create_ref_genome_from_server_location'),
    url(r'^_/ref_genomes/create_from_ncbi$',
            'main.xhr_handlers.create_ref_genome_from_ncbi'),
    url(r'^_/ref_genomes/delete$',
            'main.xhr_handlers.ref_genomes_delete'),
    url(r'^_/ref_genomes/concatenate$',
            'main.xhr_handlers.ref_genomes_concatenate'),
    url(r'^_/ref_genomes/download$',
            'main.xhr_handlers.ref_genomes_download'),
    url(r'^_/single_ref_genome$',
            'main.xhr_handlers.get_single_ref_genome'),

    url(r'^_/contigs$',
        'main.xhr_handlers.get_contigs'),
    url(r'^_/contigs/delete$',
        'main.xhr_handlers.contigs_delete'),

    ###########################################################################
    # Template XHR's
    # TODO: Replace this with client-side templating.
    ###########################################################################

    url(r'^_/templates/variant_filter_controls$',
            'main.template_xhrs.variant_filter_controls'),
    url(r'^_/templates/variant_set_list_controls$',
            'main.template_xhrs.variant_set_list_controls'),
    url(r'^_/templates/alignment_controls$',
            'main.template_xhrs.alignment_controls'),
    url(r'^_/templates/alignment_list_controls$',
            'main.template_xhrs.alignment_list_controls'),
    url(r'^_/templates/reference_genome_list_controls$',
            'main.template_xhrs.reference_genome_list_controls'),
    url(r'^_/templates/sample_list_controls$',
            'main.template_xhrs.sample_list_controls'),
    url(r'^_/templates/contig_list_controls$',
            'main.template_xhrs.contig_list_controls'),
    url(r'^_/templates/create_new_empty_variant_set$',
            'main.template_xhrs.create_new_empty_variant_set'),

    ###########################################################################
    # Jbrowse Redirect
    # For re-compiling the trackList.json before redirecting to
    # the static jbrowse/ pages.
    ###########################################################################
    url(r'^redirect_jbrowse$',
            'main.views.compile_jbrowse_and_redirect'),

)


###########################################################################
# Auth
###########################################################################

urlpatterns += patterns('',
        # The default behavior of registration is redirect to 'users/<username>'.
        # For now let's catch this request here and just redirect to '/'.
        (r'^users/', RedirectView.as_view(url='/')),
)

if settings.DEMO_MODE:
    from django.contrib.auth import views as auth_views
    urlpatterns += patterns('',
        # Automatically log users in during demo.
        url(r'^accounts/login/$',
                'main.demo_view_overrides.login_demo_account',
                name='auto_login'),
        url(r'^accounts/logout/$',
               auth_views.logout,
               {'template_name': 'registration/logout.html'},
               name='auth_logout'),
    )
else:
    urlpatterns += patterns('',
        # django-registration defaults
        # (further delgates to django.contrib.auth.url)
        url(r'^accounts/register/$',
               RegistrationViewWrapper.as_view(),
               name='registration_register'),
        url(r'^accounts/register/closed/$',
               TemplateView.as_view(template_name=
                        'registration/registration_closed.html'),
               name='registration_disallowed'),
        url(r'^accounts/', include('registration.auth_urls')),

        # The default behavior of registration is redirect to 'users/<username>'.
        # For now let's catch this request here and just redirect to '/'.
        url(r'^users/', RedirectView.as_view(url='/')),
    )



###########################################################################
# S3
###########################################################################

if settings.DEBUG:
    from django.conf.urls.static import static
    urlpatterns += static('jbrowse', document_root=settings.JBROWSE_ROOT,
            show_indexes=True)


###########################################################################
# S3
###########################################################################

if settings.S3_ENABLED:
    urlpatterns += patterns('',
        url(r'^_/projects/([\w-]+)/refgenomes/import_s3$',
                'main.xhr_handlers.import_reference_genome_s3',
                name="import_reference_genome_s3"),
        url(r'^_/projects/([\w-]+)/samples/parse_targets_file_s3$',
                'main.xhr_handlers.parse_targets_file_s3',
                name="parse_targets_file_s3"),
        url(r'^_/projects/([\w-]+)/samples/process_sample_files_s3$',
                'main.xhr_handlers.process_sample_files_s3',
                name="process_sample_files_s3"),
        url(r'^s3/signature', 'main.xhr_uploader.handle_s3',
                name="s3_signature"),
        url(r'^s3/delete', 'main.xhr_uploader.handle_s3',
                name='s3_delete'),
        url(r'^s3/success', 'main.xhr_uploader.success',
                name="s3_success")
    )

# TODO: Uncomment once we're sure we want this.
# if settings.RUNNING_ON_EC2:
#     urlpatterns += patterns('',
#         url(r'^ec2/info$', 'main.views.ec2_info_view',
#                 name="ec2_info")
#     )

###########################################################################
# Functions to run at startup
###########################################################################

from main import startup
startup.run()
