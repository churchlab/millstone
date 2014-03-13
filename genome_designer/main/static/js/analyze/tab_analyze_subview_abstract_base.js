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
  },


  /** Register control listeners. Inheriting classes should override. */
  listenToControls: function() {
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
