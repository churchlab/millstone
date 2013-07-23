/**
 * @fileoverview The home view containing search results.
 */


gd.SampleView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-samples').addClass('active');    
  },
  
  events: {
    'click #submitFormFromFile': 'handleFormSubmit',
  },
  
  handleFormSubmit: function() {
    $("#formFromFile").submit();
  },
  
});
