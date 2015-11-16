from collections import namedtuple

from django.conf import settings


Junction = namedtuple('Junction',
        ['ref', 'ref_count', 'contig', 'contig_count'])


def get_ref_jbrowse_link(contig, loc):
    sample_alignment = contig.experiment_sample_to_alignment
    # This is a workaround to using type=Dataset.TYPE.BWA_ALIGN, because this
    # file is imported by main/models.py, so importing Dataset here causes a
    # circular import error.
    # TODO(gleb): Figure out better organization
    bam_dataset = sample_alignment.dataset_set.get(
            type='BWA BAM')

    sample_bam_track = '_'.join([
            bam_dataset.internal_string(sample_alignment.experiment_sample),
            str(sample_alignment.alignment_group.uid)])

    track_labels = settings.JBROWSE_DEFAULT_TRACKS + [sample_bam_track]
    return (contig.parent_reference_genome.get_client_jbrowse_link() +
            '&loc=' + str(loc) +
            '&tracks=' + ','.join(track_labels))


def decorate_with_link_to_loc(contig, loc, text):
    return ('<a href="' + get_ref_jbrowse_link(contig, loc) +
        ',' + contig.get_contig_reads_track() +
        '" target="_blank">' + text + '</a>')


def make_html_list(li, css_class='list-unstyled'):
    return ('<ul class="' + css_class + '"><li>' +
            '</li><li>'.join(li) +
            '</li></ul>')


def create_contig_junction_links(contig, junctions):
    link_list = ['<span>&rarr;</span>'.join(
            [decorate_with_link_to_loc(contig, j.ref,
                    '%s(%s)' % (j.ref, j.ref_count)),
            '%s(%s)' % (j.contig, j.contig_count)])
            for j in map(Junction._make, junctions)]
    return make_html_list(link_list)
