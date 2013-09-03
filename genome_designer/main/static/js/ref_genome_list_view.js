/**
 * @fileoverview Reference Genome List view.
 */


gd.RefGenomeListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-refgenomes').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-ref-genome-list-view-datatable-hook'),
        objList: REF_GENOME_LIST_DATA['obj_list'],
        fieldConfig: REF_GENOME_LIST_DATA['field_config']
    });

  }
});
