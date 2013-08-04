/**
 * @fileoverview Single alignment and samples within.
 */


gd.AlignmentView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    this.decorateSidebar();

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-alignment_view-datatable-hook'),
        objList: EXPERIMENT_TO_SAMPLE_DATA['obj_list'],
        fieldConfig: EXPERIMENT_TO_SAMPLE_DATA['field_config'],
    });
    
    this.datatable.addDropdownOption('Test Option', 'test');
  },

  /** Decorate the side nav bar. */
  decorateSidebar: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    // Draw a sub-menu.
    $('#gd-sidenav-link-alignments').append(
        '<ul class="nav nav-list">' +
          '<li>... ' + this.model.get('label') + '</li>' +
        '</ul>');
  },

  events: {
    'click #gd-alignments-rerun-variants-btn': 'handleRerunVariantsClick',
  },

  /** Tells the server to re-run variant calling for this alignment group. */
  handleRerunVariantsClick: function() {
    $('#gd-alignments-rerun-variants-btn').prop('disabled', true);

    // Post to this same view for now.
    var postUrl = window.location.pathname;

    // Nothing to send as this is the only post request this page expects.
    var postData = {};

    var onSuccess = function(data) {
      $('#gd-alignments-rerun-variants-btn').text('Re-running...');
    };

    // Execute the post. Should return a redirect response.
    $.post(postUrl, JSON.stringify(postData), onSuccess, 'json');
  }
});
