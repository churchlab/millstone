/**
 * @fileoverview List of variant sets.
 */


gd.VariantSetListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-variant-sets').addClass('active');
    
    this.datatable = new gd.DataTableComponent({
        el: $('#gd-variant_set_list_view-datatable-hook'),
        objList: VARIANT_SET_LIST_DATA['obj_list'],
        fieldConfig: VARIANT_SET_LIST_DATA['field_config']
    });
  }
});
