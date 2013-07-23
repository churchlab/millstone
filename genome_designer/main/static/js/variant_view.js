/**
 * @fileoverview List of alignments.
 */


gd.VariantView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-variants').addClass('active');
  }
});
