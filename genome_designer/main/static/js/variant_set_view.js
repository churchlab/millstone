/**
 * @fileoverview Single alignment and samples within.
 */


gd.VariantSetView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    // Field config for the datatable.
    this.fieldConfig = VARIANT_TO_VARIANT_SET_DATA.field_config;

    // List of objects to display in the datatable.
    this.variantToVariantSetData = VARIANT_TO_VARIANT_SET_DATA.obj_list;

    this.render();
  },

  render: function() {
    this.decorateSidebar();

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-datatable-hook'),
        objList: this.variantToVariantSetData,
        fieldConfig: this.fieldConfig
    });

  },

  /** Decorate the side nav bar. */
  decorateSidebar: function() {
    $('#gd-sidenav-link-variant-sets').addClass('active');
  },

  events: {
    'click #gd-variant-set-view-export-as-csv': 'handleExportCsv',
  },

  handleExportCsv: function() {
    $('#gd-variant-sets-export-csv-form').submit();
  }
});
