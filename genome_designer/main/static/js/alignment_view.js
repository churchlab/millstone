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
        fieldConfig: EXPERIMENT_TO_SAMPLE_DATA['field_config']
    });
  },

  /** Decorate the side nav bar. */
  decorateSidebar: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    // Draw a sub-menu.
    $('#gd-sidenav-link-alignments').append(
        '<ul class="nav nav-list">' +
          '<li>... ' + this.model.get('label') + '</li>' +
        '</ul>');
  }
});
