/**
 * @fileoverview List of alignments.
 */


gd.GotermView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-goterms').addClass('active');
  }
});
