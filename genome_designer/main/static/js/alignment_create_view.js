/**
 * @fileoverview List of variants.
 */


gd.AlignmentCreateView = Backbone.View.extend({
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
    $('#gd-datatable-ref_genome-hook').html(
      '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered"' +
          'id="gd-ref_genome-datatable">' +
      '</table>');
    $('#gd-datatable-samples-hook').html(
      '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered"' +
          'id="gd-samples-datatable">' +
      '</table>');
    $('#gd-ref_genome-datatable').dataTable({
        'aaData': REF_GENOME_LIST_DATA['obj_list'],
        'aoColumns': REF_GENOME_LIST_DATA['field_config'],
        "bSortClasses": false,
        'sDom': "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
        'sPaginationType': 'bootstrap'
    });
    $('#gd-samples-datatable').dataTable({
        'aaData': SAMPLES_LIST_DATA['obj_list'],
        'aoColumns': SAMPLES_LIST_DATA['field_config'],
        "bSortClasses": false,
        'sPaginationType': 'bootstrap'
    });
    
  }
});
