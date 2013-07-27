/**
 * @fileoverview List of variants.
 */


gd.VariantView = Backbone.View.extend({
  el: '#gd-page-container',


  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-variants').addClass('active');
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
    $('#gd-datatable-hook').html(
        '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered"' +
            'id="gd-datatable">' +
        '</table>');
    $('#gd-datatable').dataTable({
        'aaData': VARIANT_LIST_DATA['obj_list'],
        'aoColumns': VARIANT_LIST_DATA['field_config'],
        'sDom': "<'row'<'span5'l><'span5'f><'align-right span2'C>r>t<'row'<'span6'i><'span6'p>>",
        "bSortClasses": false,
        'sPaginationType': 'bootstrap'
    });
  }
});
