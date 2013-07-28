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

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-datatable-hook'),
        objList: VARIANT_LIST_DATA['obj_list'],
        fieldConfig: VARIANT_LIST_DATA['field_config']
    });
  }
});
