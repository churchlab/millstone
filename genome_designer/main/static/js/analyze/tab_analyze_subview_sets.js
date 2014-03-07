/**
 * @fileoverview Sub-view displaying VariantSets.
 */

gd.TabAnalyzeSubviewSets = gd.TabAnalyzeSubviewAbstractBase.extend({


  /** @override */
  initialize: function() {
    // The DataTableComponent currently in view.
    this.datatableComponent = null;

    gd.TabAnalyzeSubviewAbstractBase.prototype.initialize.call(this);
  },


  /** @override */
  renderControls: function() {
    // Request the filter control html from the server for now.
    // TODO: Client-side template support.
    $.get('/_/templates/variant_set_controls', _.bind(function(response) {
      // Append the DOM.
      $('.gd-table-controls').append(response);

      // Register event listeners.
      $('#submitFormFromFile').click(
          _.bind(this.handleFormSubmitFromFile, this));
      $('#gd-variant-set-form-empty-submit').click(
          _.bind(this.handleFormSubmitEmpty, this));
    }, this));
  },


  /** Handles a click on the modal for creating a new VariantSet from file. */
  handleFormSubmitFromFile: function() {
    // Page form submit is not really the right thing here.
    alert('TODO');
  },


  /** Handles a click on the modal for creating a new empty VariantSet. */
  handleFormSubmitEmpty: function() {
    // Page form submit is not really the right thing here.
    alert('TODO');
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

    $.get('/_/sets', requestData,
        _.bind(this.handleGetVariantSetListResponse, this));
  },


  /** Handles the data fetch response. */
  handleGetVariantSetListResponse: function(response) {
    var objList = [];
    var fieldConfig = {};

    // Parse the response.
    var data = JSON.parse(response.variant_set_list_json);
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
  }
});
