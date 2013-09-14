/**
 * @fileoverview List of variants.
 */


gd.SingleVariantView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
  // Field config for the datatable.
    this.fieldConfig = MELTED_VARIANT_DATA.field_config;

    // List of objects to display in the datatable.
    this.tableData = MELTED_VARIANT_DATA.obj_list;

    this.render();
  },

  render: function() {
    this.datatable = new gd.DataTableComponent({
        el: $('#gd-datatable-hook'),
        objList: this.tableData,
        fieldConfig: this.fieldConfig
    });
  },
});
