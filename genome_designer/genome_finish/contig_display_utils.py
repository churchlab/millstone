from collections import namedtuple

Junction = namedtuple('Junction', ['ref', 'ref_count', 'contig', 'contig_count'])


def get_ref_jbrowse_link(contig, loc):
        return (contig.parent_reference_genome.get_client_jbrowse_link() +
                '&loc=' + str(loc))


def decorate_with_link_to_loc(contig, loc, text):
    return ('<a href="' + get_ref_jbrowse_link(contig, loc) +
        '&tracks=' + contig.get_contig_reads_track() +
        '" target="_blank">' + text + '</a>')


def make_html_list(li, css_class='list-unstyled'):
    return ('<ul class="' + css_class + '"><li>' +
            '</li><li>'.join(li) +
            '</li></ul>')


def create_contig_junction_links(contig, junctions):
    link_list = ['<span>&rarr;</span>'.join(
            [decorate_with_link_to_loc(contig, j.ref, '%s(%s)' % (j.ref, j.ref_count)),
            '%s(%s)' % (j.contig, j.contig_count)]) for j in map(Junction._make, junctions)]
    return make_html_list(link_list)
