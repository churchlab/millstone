/**
 * @fileoverview Sub-view displaying genes.
 */

gd.TabAnalyzeSubviewGenes = gd.TabAnalyzeSubviewAbstractBase.extend({

  /** @override */
  initialize: function() {
    // The DataTableComponent currently in view.
    this.datatableComponent = null;

    gd.TabAnalyzeSubviewAbstractBase.prototype.initialize.call(this);
  },


  /** @override */
  renderDatatable: function() {
    // Request data and load the datatable.
    this.createDatatableComponent();
  },


  /** Redraws the datatable based on the selection. */
  createDatatableComponent: function() {
    var requestData = {
      'refGenomeUid': this.model.get('refGenomeUid'),
    };

    $.get('/_/genes', requestData,
        _.bind(this.handleGetGenesResponse, this));
  },


  /** Handles the data fetch response. */
  handleGetGenesResponse: function(response) {
    var objList = [];
    var fieldConfig = {};

    // Parse the response.
    var data = JSON.parse(response.geneList);
    if (data.obj_list.length) {
      objList = data.obj_list;
      fieldConfig = data.field_config;
    }

    // Draw the datatable.
    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-datatable-hook'),
        objList: objList,
        fieldConfig: fieldConfig
    });
  },
});
