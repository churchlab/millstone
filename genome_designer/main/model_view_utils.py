"""
model_view_utils - Functions to decorate model_view output before sending
to datatables. 
"""

from collections import defaultdict
from math import floor
from itertools import chain

from django.core.urlresolvers import reverse

import settings
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__CHROMOSOME
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__ES_UID
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__POSITION


def get_jbrowse_track_names(reference_genome):
    """
    This is a master list of all jbrowse track names for a reference genome,
    ordered by type. 

    TODO: We should generate this programmatically, maybe store it in the db.
    """
    jbrowse_track_names = {}
    jbrowse_track_names['vcf'] = []

    alignment_groups = reference_genome.alignmentgroup_set.all()

    # VCF
    for ag in alignment_groups:      
        jbrowse_track_names['vcf'].append('_'.join([
            str(ag.uid),
            'SNPEFF_VCF']))

    # BAM
    jbrowse_track_names['bam'] = ['_%s_%s' % (
            'BWA_BAM', ag.uid) for ag in alignment_groups]

    # BAM_COVERAGE
    jbrowse_track_names['bam_coverage'] = ['_%s_%s' % (
            'BWA_BAM_COVERAGE', ag.uid) for ag in alignment_groups]

    # CALLABLE_LOCI_BED
    jbrowse_track_names['callable_loci_bed'] = ['_%s_%s' % (
            'FLAGGED_REGIONS_BED', ag.uid) for ag in alignment_groups]

    return jbrowse_track_names


# HELPER FXNS FOR GENERATING MODEL VIEW LINKS ============================

def create_single_variant_page_link_for_variant_object(variant_as_dict,
        reference_genome, alignment_group):
    """Constructs the label as an anchor that links to the single variant view.
    """
    # Generate link to Analyze view for this single variant.
    root_href = reverse('main.views.tab_root_analyze',
            args=(reference_genome.project.uid,
                    alignment_group.uid, 'variants'))
    filter_part = '?filter=UID=%s&melt=1' % (variant_as_dict['UID'],)
    full_href = root_href + filter_part

    return full_href


def create_jbrowse_link_for_variant_object(variant_as_dict, reference_genome,
        track_strings):
    """
    Constructs a JBrowse link for the Variant. Adds tracks passed in; 
    DNA and genbank annotation tracks are in by default. 
    """
    assert MELTED_SCHEMA_KEY__POSITION in variant_as_dict
    position = variant_as_dict[MELTED_SCHEMA_KEY__POSITION]
    ref_genome_jbrowse_link = reference_genome.get_client_jbrowse_link()
    
    location_str = '..'.join([str(i) for i in [
            position - int(floor(settings.JBROWSE_DEFAULT_VIEW_WINDOW/2)),
            position + int(floor(settings.JBROWSE_DEFAULT_VIEW_WINDOW/2))]])

    location_param = '&loc=' + variant_as_dict[MELTED_SCHEMA_KEY__CHROMOSOME] + ':' + str(location_str)

    tracks = [] + settings.JBROWSE_DEFAULT_TRACKS

    # add all specified tracks to this view
    tracks.extend(track_strings)

    tracks_param = '&tracks=' + ','.join(tracks)

    full_href = ref_genome_jbrowse_link + location_param + tracks_param

    return full_href

def create_variant_links_field(variant_as_dict, reference_genome,
        alignment_group, jbrowse_track_names):
    """
    Create a list of icon links for the variant datatable view.
    """

    # list of all experiment samples for this variant row
    es_list = []

    # add all experiment samples to get their track strings below
    if MELTED_SCHEMA_KEY__ES_UID in variant_as_dict:
        es_field = variant_as_dict[MELTED_SCHEMA_KEY__ES_UID]

        if es_field is None: es_list = []

        # melted view, one experiment sample:
        elif isinstance(es_field, basestring) and es_field is not None:
            es_list.append(es_field)

        else:
            es_list = [es for es in es_field if es is not None]

    # BAM JBROWSE
    jbrowse_bam_tracks = list(chain.from_iterable([
            jbrowse_track_names['vcf'] + 
            [es + s for s in jbrowse_track_names['bam']] + 
            [es + s for s in jbrowse_track_names['callable_loci_bed']]
                    for es in es_list]))

    jbrowse_bam_href = create_jbrowse_link_for_variant_object(
            variant_as_dict, reference_genome, jbrowse_bam_tracks)

    # BAM COVERAGE JBROWSE
    jbrowse_bam_coverage_tracks = list(chain.from_iterable([
            jbrowse_track_names['vcf'] + 
            [es + s for s in jbrowse_track_names['bam_coverage']] + 
            [es + s for s in jbrowse_track_names['callable_loci_bed']]
                    for es in es_list]))

    jbrowse_bam_coverage_href = create_jbrowse_link_for_variant_object(
            variant_as_dict, reference_genome, jbrowse_bam_coverage_tracks)

    single_variant_view_href = create_single_variant_page_link_for_variant_object(
            variant_as_dict, reference_genome,
            alignment_group)

    buttons = [{
            'href': single_variant_view_href,
            'glyph': 'glyphicon-search',
            'title': 'single variant view',
            'target': '_self'
        },
        { 
            'href': jbrowse_bam_href, 
            'glyph': 'glyphicon-sort-by-attributes',
            'title': 'BAM alignment',
            'max': settings.JBROWSE_MAX_ALIGN_TRACKS
        },
        {
            'href': jbrowse_bam_coverage_href,
            'glyph': 'glyphicon-stats',
            'title': 'BAM coverage',
            'max': settings.JBROWSE_MAX_COVERAGE_TRACKS
        }]

    all_buttons_html = []
    for button in buttons:

        # If there are too many samples, don't display the tracks, 
        # and gray-out the icon
        if 'max' in button and len(es_list) > button['max']:
            button['href'] = create_jbrowse_link_for_variant_object(
                    variant_as_dict, 
                    reference_genome, 
                    jbrowse_track_names['vcf'])
            button['title'] += ' (too many samples)'
            button['glyph'] += ' disabled'

        button_html = ('<a target=' + button.get('target','"_blank"') +
            ' href="' + button.get('href','#') + '"' +
            ' title="' + button.get('title','') + '">' +
            '<span class="glyphicon ' + button['glyph'] + 
            '"></span></a>')

        all_buttons_html.append(button_html)

    return ' '.join(all_buttons_html)

