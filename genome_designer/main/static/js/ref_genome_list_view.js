/**
 * @fileoverview Reference Genome List view.
 */


gd.RefGenomeListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-refgenomes').addClass('active');

    this.redrawDatatable();
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    var requestData = {
        projectUid: this.model.get('uid'),
    };

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-ref-genome-list-view-datatable-hook'),
        serverTarget: '/_/ref_genomes',
        controlsTemplate: '/_/templates/reference_genome_list_controls',
        requestData: requestData,
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  decorateControls: function() {
    if (this.refGenomeControlsComponent) {
      this.refGenomeControlsComponent.destroy();
    }

    this.refGenomeControlsComponent = new gd.RefGenomeControlsComponent({
      el: '#gd-ref-genome-list-view-datatable-hook-control',
      datatableComponent: this.datatableComponent
    });

    this.listenTo(this.refGenomeControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  }
});
