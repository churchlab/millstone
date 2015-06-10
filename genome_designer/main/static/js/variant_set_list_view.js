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
    this.redrawDatatable();
  },

  decorateControls: function() {
    this.variantsetsControlsComponent = new gd.VariantsetsControlsComponent({
      el: '#gd-variant_set_list_view-datatable-hook-control',
      model: this.model,
      datatableComponent: this.datatableComponent
    });

    this.listenTo(this.variantsetsControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    if (this.variantsetsControlsComponent) {
      this.variantsetsControlsComponent.destroy();
    }

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-variant_set_list_view-datatable-hook'),
        serverTarget: '/_/sets',
        controlsTemplate: '/_/templates/variant_set_list_controls',
        requestData: {projectUid: this.model.get('uid')}
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  }  
});
