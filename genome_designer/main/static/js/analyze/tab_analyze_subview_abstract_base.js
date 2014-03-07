/**
 * @fileoverview Abstract base class for subviews of the Analyze tab.
 */

gd.TabAnalyzeSubviewAbstractBase = Backbone.View.extend({

  /** Override. */
  initialize: function() {
    this.render();
  },


  /** Override. */
  render: function() {
    this.renderDatatable();
    // This needs to happen afterwards since they are now inside the table. 
    this.renderControls();
  },


  /** Render the controls. Inheriting classes should override. */
  renderControls: function() {
    $('.gd-table-controls').append(
        '<h3>Controls go here</h3>');
  },


  /** Render the DataTable. Inheriting classes should override. */
  renderDatatable: function() {
    $('#gd-datatable-hook').append(
        '<h3>Datatable goes here</h3>');
  },


  /** Properly clean up this component. Inheriting classes should override. */
  destroy: function() {
  }
});
