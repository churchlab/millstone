/**
 * @fileoverview List of alignments.
 */


gd.AlignmentListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-alignment_list_view-datatable-hook'),
        serverTarget: '/_/alignmentgroups',
        controlsTemplate: '/_/templates/alignment_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });
  }
});
