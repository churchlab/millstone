/**
 * @fileoverview List of alignments.
 */


gd.AlignmentListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
    this.decorate_new_button();
  },

  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-alignment_list_view-datatable-hook'),
        objList: ALIGNMENT_LIST_DATA['obj_list'],
        fieldConfig: ALIGNMENT_LIST_DATA['field_config']
    });
  },

  decorate_new_button: function() {

    $("div.gd-new-button").html(
      '<a href=' + NEW_ALIGNMENT_LINK + '>' +
        '<button type="submit" class="btn btn-primary">New...</button>' +
      '</a>'
    );
  }

});
