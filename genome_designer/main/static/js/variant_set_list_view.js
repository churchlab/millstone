/**
 * @fileoverview List of variant sets.
 */


gd.VariantSetListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
    this.decorate_new_button();
  },

  render: function() {
    $('#gd-sidenav-link-variant-sets').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-variant_set_list_view-datatable-hook'),
        objList: VARIANT_SET_LIST_DATA['obj_list'],
        fieldConfig: VARIANT_SET_LIST_DATA['field_config']
    });
  },

  events: {
    'click #submitFormFromFile': 'handleFormSubmitFromFile',
    'click #gd-variant-set-form-empty-submit': 'handleFormSubmitEmpty'
  },

  handleFormSubmitFromFile: function() {
    $("#formFromFile").submit();
  },

  handleFormSubmitEmpty: function() {
    $("#gd-variant-set-form-empty").submit();
  },

  decorate_new_button: function() {

    $("div.gd-new-button").html(
      '<div class="btn-group">' +
      '  <a class="btn dropdown-toggle btn-primary" data-toggle="dropdown" href="#">' +
      '    Create New' +
      '    <span class="caret"></span>' +
      '  </a>' +
      '  <ul class="dropdown-menu">' +
      '    <li role="presentation"><a role="menuitem" tabindex="-1" href="#modalFromFile" data-toggle="modal">From File...</a></li>' +
      '    <li role="presentation"><a role="menuitem" tabindex="-1" href="#modalEmpty" data-toggle="modal">Empty</a></li>' +
      '  </ul>' +
      '</div>'
    );
  }
});
