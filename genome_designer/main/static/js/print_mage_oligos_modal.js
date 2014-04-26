/**
 * Modal that allows choosing options for printing MAGE oligos.
 */


gd.PrintMageOligosModal = Backbone.View.extend({

  initialize: function() {
    this.render();
  },

  render: function() {
    this.listenToControls();
  },

  listenToControls: function() {
    $('#gd-variant-sets-print-mage-oligos-submit').click(
        _.bind(this.handlePrintMageOligos, this));
  },

  handlePrintMageOligos: function() {
    $('#gd-variant-set-print-mage-oligos-form').submit();
  }
});
