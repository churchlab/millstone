/**
 * @fileoverview Single alignment and samples within.
 */


gd.VariantSetView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    this.decorateSidebar();

    var filterString = 'variant_set_uid = ' + this.model.get('uid');
    this.variantsTableComponent = new gd.VariantsTableComponent({
      model: this.model,
      filterString: filterString,
      filterDisabled: true,
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
