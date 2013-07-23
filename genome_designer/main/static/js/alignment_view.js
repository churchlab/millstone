/**
 * @fileoverview List of alignments.
 */


gd.AlignmentView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');
  }
});
