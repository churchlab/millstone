/**
 * @fileoverview Single alignment and samples within.
 */


gd.VariantSetView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    this.decorateSidebar();

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-variant_set_view-datatable-hook'),
        objList: VARIANT_TO_VARIANT_SET_DATA['obj_list'],
        fieldConfig: VARIANT_TO_VARIANT_SET_DATA['field_config']
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
