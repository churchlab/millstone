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
    this.updateDatatable();
  },

  /**
   * Updates the datatable view based on the data available.
   *
   * @param {array} data List of objects to display
   * @param {array} column_config List of column config objects. These must
   *    have the keys (NOTE camel-case):
   *        mData: key corresponding to key in data.
   *        sTitle: title for the column.
   */
  updateDatatable: function(data, column_config) {
    $('#gd-alignment_view-datatable-hook').html(
        '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered"' +
            'id="gd-alignment_view-datatable">' +
        '</table>');
    $('#gd-alignment_view-datatable').dataTable({
        'aaData': ALIGNMENT_LIST_DATA['obj_list'],
        'aoColumns': ALIGNMENT_LIST_DATA['field_config'],
        "bSortClasses": false,
        'sPaginationType': 'bootstrap'
    });
  }

});