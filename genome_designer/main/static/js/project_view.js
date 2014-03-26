/**
 * @fileoverview View of a specific project.
 */


gd.ProjectView = Backbone.View.extend({
  el: '#gd-page-container',

  events: {
    'click #gd-projects-delete-btn': 'handleDeleteClick',
    'click #gd-projects-export-btn': 'handleExportClick'
  },

  handleDeleteClick: function() {
    var agree = confirm("Are you sure you want to delete ths project?");
    if (agree) {
      var postUrl = '/projects/' + this.model.get('uid') + '/delete'

      /** Callback upon successful delete. */
      var onSuccess = function(data) {
        // Redirect according to the server response.
        window.location.href = data.redirect;
      };

      $.post(postUrl, onSuccess, 'json');
    }
  },

  handleExportClick: function() {
    alert('Coming soon.');
  }
});
