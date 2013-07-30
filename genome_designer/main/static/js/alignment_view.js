/**
 * @fileoverview List of alignments.
 */


gd.AlignmentView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-alignment_view-datatable-hook'),
        objList: ALIGNMENT_LIST_DATA['obj_list'],
        fieldConfig: ALIGNMENT_LIST_DATA['field_config']
    });
  },
});
