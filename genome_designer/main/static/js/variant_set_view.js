/**
 * @fileoverview List of alignments.
 */


gd.VariantSetView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-variant-sets').addClass('active');
  }
});
