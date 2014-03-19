/**
 * @fileoverview Single ReferenceGenome view.
 */


gd.RefGenomeView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-refgenomes').addClass('active');
  },
});
