/**
 * @fileoverview Sub-view displaying Contigs.
 */


gd.TabAnalyzeSubviewContigs = gd.TabAnalyzeSubviewAbstractBase.extend({

  /** @override */
  initialize: function() {
    this.render();
  },

  render: function() {
    this.redrawDatatable();
    this.redrawSamplesDatatable();
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    var requestData = {
        refGenomeUid: this.model.attributes.refGenomeUid,
        alignmentGroupUid: this.model.attributes.alignmentGroupUid
    };

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-datatable-hook-contigs'),
        serverTarget: '/_/contigs',
        controlsTemplate: '/_/templates/contig_list_controls',
        requestData: requestData,
    });

    this.listenToOnce(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  redrawSamplesDatatable: function() {
    if (this.samplesDatatableComponent) {
      this.samplesDatatableComponent.destroy();
    }

    var samplesRequestData = {
        alignmentGroupUid: this.model.attributes.alignmentGroupUid
    };

    this.samplesDatatableComponent = new gd.DataTableComponent({
        el: $('#samples-datatable-hook'),
        serverTarget: '/_/sample_alignments_for_assembly',
        controlsTemplate: '/_/templates/sample_alignment_list_for_assembly_controls',
        requestData: samplesRequestData,
        extraDatatableParams: {'iDisplayLength': 10}
    });

    this.listenToOnce(this.samplesDatatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.samplesDecorateControls, this));
  },

  decorateControls: function() {
    if (this.contigControlsComponent) {
      this.contigControlsComponent.destroy();
    }

    this.contigControlsComponent = new gd.ContigControlsComponent({
      el: '#gd-datatable-hook-control',
      datatableComponent: this.datatableComponent,
      alignmentGroupUid: this.model.attributes.alignmentGroupUid
    });

    this.listenToOnce(this.contigControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  },

  samplesDecorateControls: function() {

    this.sampleAssemblyControlsComponent = new gd.SampleAssemblyControlsComponent({
      el: '#samples-datatable-hook-control',
      datatableComponent: this.samplesDatatableComponent,
      alignmentGroupUid: this.model.attributes.alignmentGroupUid
    });

    this.listenToOnce(this.sampleAssemblyControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawSamplesDatatable, this));
  }
});
