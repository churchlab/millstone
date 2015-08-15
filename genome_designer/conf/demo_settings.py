"""
Settings for DEMO_MODE.

Must set DEMO_MODE = True in local_settings.py.
"""

# Views that are visible in demo mode.
DEMO_SAFE_VIEWS = [
    'django.contrib.auth.views.logout',
    'main.demo_view_overrides.login_demo_account',
    'main.template_xhrs.alignment_controls',
    'main.template_xhrs.reference_genome_list_controls',
    'main.template_xhrs.sample_list_controls',
    'main.template_xhrs.variant_filter_controls',
    'main.template_xhrs.variant_set_list_controls',
    'main.views.alignment_list_view',
    'main.views.alignment_view',
    'main.views.compile_jbrowse_and_redirect',
    'main.views.home_view',
    'main.views.project_list_view',
    'main.views.project_view',
    'main.views.reference_genome_list_view',
    'main.views.reference_genome_view',
    'main.views.sample_alignment_error_view',
    'main.views.sample_list_view',
    'main.views.tab_root_analyze',
    'main.views.variant_set_list_view',
    'main.views.variant_set_view',
    'main.xhr_handlers.get_alignment_groups',
    'main.xhr_handlers.get_gene_list',
    'main.xhr_handlers.get_ref_genomes',
    'main.xhr_handlers.get_samples',
    'main.xhr_handlers.get_variant_list',
    'main.xhr_handlers.get_variant_set_list',
    'main.xhr_handlers.is_materialized_view_valid',
    'main.xhr_handlers.refresh_materialized_variant_table',
]
