/**
 * @fileoverview Reference Genome List view.
 */


gd.GenomeFinishView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-genome-finish').addClass('active');

    this.redrawDatatable();
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-genome-finish-view-datatable-hook'),
        serverTarget: '/_/ref_genomes',
        controlsTemplate: '/_/templates/reference_genome_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  decorateControls: function() {
    this.refGenomeControlsComponent = new gd.RefGenomeControlsComponent({
      el: '#gd-genome-finish-view-datatable-hook-control',
      datatableComponent: this.datatableComponent
    });

    this.listenTo(this.refGenomeControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  }
});
