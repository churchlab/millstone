/**
 * @fileoverview List of variants.
 */


gd.VariantListView = Backbone.View.extend({
  el: '#gd-page-container',


  initialize: function() {
    this.render();
  },


  render: function() {
    $('#gd-sidenav-link-variants').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-datatable-hook'),
        objList: VARIANT_LIST_DATA['obj_list'],
        fieldConfig: VARIANT_LIST_DATA['field_config']
    });

    this.draw_dropdowns();

  },

  draw_dropdowns: function() {
    // Create 'add selected to set' button, display list of available sets
    this.addVariantSetDropdownSubmenu(
      VARIANT_SET_LIST_DATA['obj_list'], 'add', 'Add to set');
    this.addVariantSetDropdownSubmenu(
      VARIANT_SET_LIST_DATA['obj_list'], 'remove', 'Remove from set');

  },

  events: {
    'click .gd-variant-set-action': 'handleVariantSetClick',
  },

  /*
   * Create a dropdown submenu with a list of all variant sets, to which
   * selected variants will be added.
   */
  addVariantSetDropdownSubmenu: function(allVariantSets, varSetAction,
      varSetActionText) {

    var setListHTML = (
      '<li class="dropdown-submenu pull-left">' +
        '<a tabindex="-1" href="#">' +
          varSetActionText +
        '</a>' +
        '<ul class="dropdown-menu">' +
          _.map(allVariantSets,
            function(variantSet) {
              return this.addIndividualSetListItem(variantSet, varSetAction);
            }, this).join('') +
        '</ul>' +
      '</li>');

    this.datatable.addDropdownOption(setListHTML, '');
  },

  /*
   * Create a single html list object for a sub-dropdown set option.
   *
   * @param variantSet  the variant set dictionary object from the
   *  django adapter
   * @param variantSetAction  the action to perform, either 'add' or 'remove'
   */
  addIndividualSetListItem: function(variantSet, variantSetAction) {
    return (
      '<li><a tabindex="-1" class="gd-variant-set-action" ' +
        'data-variant-set-uid="' + variantSet.uid + '" ' +
        'data-variant-set-action="' + variantSetAction +
        '" href="#">' + variantSet.label + '</a></li>');
  },

  handleVariantSetClick: function(ev) {

    // Post to this same view for now.
    var postUrl = window.location.pathname;

    // Grab the selected rows.
    var postData = {
        variantUidList: this.datatable.getCheckedRowUids(),
        variantSetAction: $(ev.target).data('variant-set-action'),
        variantSetUid: $(ev.target).data('variant-set-uid')
    };

    var validationResult = this.validatePostData(postData);
    if (!validationResult.is_success) {
      $('#gd-variant-set-action-submit-alert-msg').text(validationResult.error_msg);
      $('#gd-variant-set-action-submit-alert').removeClass('alert-warn alert-info alert-error');
      $('#gd-variant-set-action-submit-alert').addClass('alert-error');
      $('#gd-variant-set-action-submit-alert').show();
      return;
    }

    var onSuccess = _.bind(function(response) {
      // Repopulate json data from response
      VARIANT_LIST_DATA = JSON.parse(response.variant_list_json);
      VARIANT_SET_LIST_DATA = JSON.parse(response.variant_set_list_json);

      // Draw new message from response
      $('#gd-variant-set-action-submit-alert-msg').text(response.alert_msg);
      $('#gd-variant-set-action-submit-alert').removeClass('alert-warn alert-info alert-error');
      $('#gd-variant-set-action-submit-alert').addClass('alert-' +
        response.alert_type);
      $('#gd-variant-set-action-submit-alert').show();

      this.datatable.update(VARIANT_LIST_DATA.obj_list);
      this.draw_dropdowns();
    }, this);

    // Execute the post. Should return a redirect response.
    $.post(postUrl, JSON.stringify(postData), onSuccess);
  },

  /**
   * Validate the post data before submitting.
   * @param {object} postData The creation data to be posted to the server.
   * @return {object} Response with keys:
   *     * is_success {boolean} Whether validation succeded.
   *     * error_msg {string} Human-readable description of the error.
   */
  validatePostData: function(postData) {
    if (!postData.variantUidList.length) {
      return {
          is_success: false,
          error_msg: 'Please select at least one variant.'
      };
    }

    return {
          is_success: true,
          error_msg: ''
    };
  }
});
