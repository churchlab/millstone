/**
 * @fileoverview List of alignments.
 */


gd.GeneView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-genes').addClass('active');
  }
});
