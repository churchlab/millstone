"""
model_view_utils - Functions to decorate model_view output before sending
to datatables.
"""

from itertools import chain
from itertools import groupby
import os
import re

from django.conf import settings
from django.core.urlresolvers import reverse

from main.models import Variant
from main.models import VariantAlternate
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__ALT
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__CHROMOSOME
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__ES_UID
from variants.melted_variant_schema import MELTED_SCHEMA_KEY__HET
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
        reference_genome, alignment_group=None):
    """Constructs the label as an anchor that links to the single variant view.
    """
    reverse_args = [reference_genome.project.uid]
    if alignment_group is not None:
        reverse_args += [alignment_group.uid, 'variants']
    root_href = reverse('main.views.tab_root_analyze', args=reverse_args)
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

    # HACK(gleb): JBrowse parses the reference .fasta file to use the first
    # word as the name of the reference genome, while we store the entire
    # string in the chromosome name. For now, imitate this.
    ref_genome_name = variant_as_dict[MELTED_SCHEMA_KEY__CHROMOSOME].split()[0]

    # Create the location string.
    window_start = position - settings.JBROWSE_DEFAULT_VIEW_WINDOW / 2
    if (variant_as_dict.get('VA_DATA', None) is not None and
            'INFO_SVLEN' in variant_as_dict['VA_DATA']):
        svlen = abs(int(variant_as_dict['VA_DATA']['INFO_SVLEN']))
        window_end = position + svlen + settings.JBROWSE_DEFAULT_VIEW_WINDOW / 2
    else:
        window_end = position + settings.JBROWSE_DEFAULT_VIEW_WINDOW / 2
    location_str = '{start}..{end}'.format(start=window_start, end=window_end)

    location_param = '&loc={ref}:{loc_str}'.format(
            ref=ref_genome_name,
            loc_str=location_str)

    tracks = [] + settings.JBROWSE_DEFAULT_TRACKS

    # add all specified tracks to this view
    tracks.extend(track_strings)

    tracks_param = '&tracks=' + ','.join(tracks)

    full_href = ref_genome_jbrowse_link + location_param + tracks_param

    return full_href


def create_variant_links_field(variant_as_dict, reference_genome,
        jbrowse_track_names, alignment_group=None):
    """
    Create a list of icon links for the variant datatable view.
    """

    # list of all experiment samples for this variant row
    es_list = []

    # add all experiment samples to get their track strings below
    if MELTED_SCHEMA_KEY__ES_UID in variant_as_dict:
        es_field = variant_as_dict[MELTED_SCHEMA_KEY__ES_UID]

        if es_field is None:
            es_list = []

        # melted view, one experiment sample:
        elif isinstance(es_field, basestring) and es_field is not None:
            es_list.append(es_field)

        else:
            es_list = [es for es in es_field if es is not None]

    variant = Variant.objects.get(uid=variant_as_dict['UID'])
    variant_specific_tracks = variant.variant_specific_tracks

    # BAM JBROWSE
    jbrowse_bam_tracks = list(chain.from_iterable([
            jbrowse_track_names['vcf'] +
            variant_specific_tracks['alignment'] +
            [es + s for s in jbrowse_track_names['bam']] +
            [es + s for s in jbrowse_track_names['callable_loci_bed']]
                    for es in es_list]))

    jbrowse_bam_href = create_jbrowse_link_for_variant_object(
            variant_as_dict, reference_genome, jbrowse_bam_tracks)

    # BAM COVERAGE JBROWSE
    jbrowse_bam_coverage_tracks = list(chain.from_iterable([
            jbrowse_track_names['vcf'] +
            variant_specific_tracks['coverage'] +
            [es + s for s in jbrowse_track_names['bam_coverage']] +
            [es + s for s in jbrowse_track_names['callable_loci_bed']]
                    for es in es_list]))

    jbrowse_bam_coverage_href = create_jbrowse_link_for_variant_object(
            variant_as_dict, reference_genome, jbrowse_bam_coverage_tracks)

    single_variant_view_href = \
            create_single_variant_page_link_for_variant_object(
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

        button_html = ('<a target=' + button.get('target', '"_blank"') +
            ' href="' + button.get('href', '#') + '"' +
            ' title="' + button.get('title', '') + '">' +
            '<span class="glyphicon ' + button['glyph'] +
            '"></span></a>')

        all_buttons_html.append(button_html)

    return ' '.join(all_buttons_html)


def create_alt_flag_field(variant_as_dict, melted, maybe_dec):
    """Display a small badge if a variant is het.
    """
    marginal_set_classes = 'gd-warn-set-badge outline'

    ve_data = variant_as_dict['VE_DATA']

    # We have two different code paths depending on whether we are rendering
    # cast or melted view. Eventualy, we should do a better job consolidating.
    if not melted:
        hets = []
        for var in ve_data:
            try:
                hets.append(var.get(MELTED_SCHEMA_KEY__HET, False))
            except:
                hets.append(False)

        processed_alts = sorted(filter(lambda alt_het: alt_het[0],
                zip(variant_as_dict[MELTED_SCHEMA_KEY__ALT], hets)))
        alt_counts = groupby(
                processed_alts, lambda alt_het: alt_het[0])

        alt_strs = []
        for alt, group in alt_counts:

            # Handle SVs of type e.g. <DEL> or <INV>. These are rendered funny
            # in html so just drop brackets here rather than mucking around
            # with escaping html.
            maybe_surrounding_brackets_match = re.match(r'<([\w]+)>', alt)
            if maybe_surrounding_brackets_match:
                alt = maybe_surrounding_brackets_match.group(1)

            # Check if BND type and add wrapper of the form BND(<...>)
            if re.match('N]', alt) or re.search('\[N', alt):
                alt = 'BND({orig_alt})'.format(orig_alt=alt)

            # Special handling in case of long alt.
            alt = maybe_handle_long_alt(alt)

            group = list(group)
            num_het = sum([alt_het[1] for alt_het in group])
            alt_string = ' %s (%d)' % (alt, len(list(group)) - maybe_dec)
            if num_het:
                alt_string += (
                        ' <span class="%s" ' +
                        'title="%d Marginal calls (IS_HET=TRUE)">' +
                        '&frac12;</span>') % (
                                marginal_set_classes,
                                num_het)

            alt_strs.append(alt_string)

        value = ' | '.join(alt_strs)

    else:
        value = variant_as_dict[MELTED_SCHEMA_KEY__ALT]
        if value is None:
            return None

        maybe_surrounding_brackets_match = re.match(r'<([\w]+)>', value)
        if maybe_surrounding_brackets_match:
            value = maybe_surrounding_brackets_match.group(1)

        # Special handling in case of long alt.
        value = maybe_handle_long_alt(value)

        if ve_data and ve_data.get(MELTED_SCHEMA_KEY__HET, False):
            value += (
                    ' <span class="%s" ' +
                    'title="Marginal call (IS_HET=TRUE)">' +
                    '&frac12;</span>') % marginal_set_classes
    return value


def maybe_handle_long_alt(alt):
    """If long alt, wrap in href appropriately.
    """
    # If this is a long alt, it will have a hash of the alt which
    # corresponds to the file on disk in which the full value is located.
    # We create an href that when clicked will result in fetching the
    # full alt from disk, guided by data-alt-hash html attribute.
    # NOTE: To be backwards compatible, we still allow long alts that
    # are stored in the database, but don't send them ot the frontend.
    # Instead we just serve up an href where data-alt-hash is empty.
    maybe_long_alt_match = VariantAlternate.LONG_ALT_REGEX.match(alt)
    if maybe_long_alt_match or len(alt) > 10:
        if maybe_long_alt_match:
            alt_value = 'LONG:{size}bp'.format(
                    size=maybe_long_alt_match.group('size'))
            alt_hash = maybe_long_alt_match.group('hash')
        else:
            alt_value = 'LONG:{size}bp'.format(size=len(alt))
            alt_hash = ''

        alt = (
                '<a '
                'class="gd-id-variants-long-alt" '
                'data-alt-hash="{alt_hash}" href="#">{alt_value}</a>'.format(
                        alt_hash=alt_hash,
                        alt_value=alt_value))
    return alt
