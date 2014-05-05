
/**
 * @fileoverview Abstract base class that provides common methods used by
 *     controls components.
 */


/** Base component with common functions for datatable controls. */
gd.DataTableControlsComponent = Backbone.View.extend({
  initialize: function() {
    // Handle to the associated DataTableComponent.
    this.datatableComponent = this.options.datatableComponent;
  },

  /**
   * Add a dropdown option to the datatable.
   */
  addDropdownOption: function (html) {
    var rendered = '<li><a href="#">' + html + '</a></li>';
    $('#' + this.datatableComponent.datatableId + '-dropdown').append(rendered);
  },
})
