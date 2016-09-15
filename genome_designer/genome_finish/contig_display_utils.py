from utils import namedtuple_with_defaults

from django.conf import settings


Junction = namedtuple_with_defaults('Junction',
        ['ref', 'ref_count', 'contig', 'contig_count', 'annotation'],
        [None, -1, None, -1, []])


def get_ref_jbrowse_link(contig, loc):
    return '#'

    # TODO(dbg): Maybe fix.
    # sample_alignment = contig.experiment_sample_to_alignment
    # # This is a workaround to using type=Dataset.TYPE.BWA_ALIGN, because this
    # # file is imported by main/models.py, so importing Dataset here causes a
    # # circular import error.
    # # TODO(gleb): Figure out better organization
    # bam_dataset = sample_alignment.dataset_set.get(
    #         type='BWA BAM')

    # sample_bam_track = '_'.join([
    #         bam_dataset.internal_string(sample_alignment.experiment_sample),
    #         str(sample_alignment.alignment_group.uid)])

    # track_labels = settings.JBROWSE_DEFAULT_TRACKS + [sample_bam_track]
    # return (contig.parent_reference_genome.get_client_jbrowse_link() +
    #         '&loc=' + str(loc) +
    #         '&tracks=' + ','.join(track_labels))


def decorate_with_link_to_loc(contig, loc, text):
    return ('<a href="' + get_ref_jbrowse_link(contig, loc) +
        # TODO(dbgoodman): Add back once JBrowse contig features re-enabled.
        # ',' + contig.get_contig_reads_track() +
        # '" target="_blank">' + text + '</a>')
        '">' + text + '</a>')


def make_html_list(li, css_class='list-unstyled'):
    return ('<ul class="' + css_class + '"><li>' +
            '</li><li>'.join(li) +
            '</li></ul>')


def create_contig_junction_links(contig, junctions):
    link_list = []

    junction_obj_list = [Junction._make(j) for j in junctions]
    for j in junction_obj_list:

        annotation_names = []
        for feature in j.annotation:

            # remove 'insertion sequence'
            if feature.startswith('insertion sequence:'):
                feature = feature[len('insertion sequence:'):]

            # skip '<unknown>'
            if '<unknown>' in feature:
                continue

            annotation_names.append(feature)

        if annotation_names:
            annotation_string = '({})'.format(
                ', '.join(annotation_names))
            ref_text = ' '.join([annotation_string, str(j.ref)])
        else:
            ref_text = str(j.ref)

        ref_text_w_link = decorate_with_link_to_loc(
                contig= contig,
                loc= j.ref,
                text= ref_text)

        contig_text = '%s' % (j.contig)

        link_list.append('<span>&rarr;</span>'.join(
                [ref_text_w_link, contig_text]))

    return make_html_list(link_list)
