/**
 * @fileoverview Base object for views related to analysis.
 */

gd.TabAnalyzeBaseView = Backbone.View.extend(
/** Prototype properties */
{
  /** Element that Backbone attches this component to. */
  el: '#gd-page-container',


  /** Backbone initialize method. */
  initialize: function() {
    // The DataTableComponent currently in view.
    this.datatableComponent = null;

    // Current list of variants being displayed.
    this.variantList = null;

    // Datatable field config.
    this.fieldConfig = null;

    // List of variant sets belonging to the user.
    this.variantSetList = null;

    // Map of variant keys.
    this.variantKeyMap = null;

    // Keys that should be visible.
    // When null, we are just showing the default.
    this.visibleKeyNames = null;

    // Perform rendering and registering listeners.
    this.render();
  },

  render: function() {
    // Manually listen to change events.
    $('#gd-analyze-select-ref-genome').change(
        _.bind(this.handleRefGenomeSelect, this));
  },


  /** Backbone sugar for registering event listeners. */
  events: {
    'click .gd-variant-set-action': 'handleVariantSetClick',
    'click #gd-filter-box-apply-btn': 'handleApplyFilterClick',
    'click #gd-filter-field-select-btn': 'handleShowFieldSelect'
  },


  /**
   * Handles a Reference Genome being selected from the dropdown menu.
   *
   * This method does an initial fetch of the data, and (re)-draws the
   * DataTable component. Subsequent requests given the same filter are
   * handled through a pagination dance between the DataTables component
   * and server-side support.
   *
   * @param {Event} e Event object if this is being called in response to an
   *     event. Ignored if opt_refGenomeUid is provided.
   * @param {string=} opt_refGenomeUid Optional ref genome string explicitly
   *     provided.  Useful for debug.
   */
  handleRefGenomeSelect: function(e, opt_refGenomeUid) {
    var refGenomeUid = opt_refGenomeUid || $(e.target).val();
    this.model.set('refGenomeUid', refGenomeUid);

    // Show the entity select (i.e. Variants, VariantSets, etc.)
    // TODO: Implement support for switching among entities.
    $('#gd-analyze-select-search-entity').fadeIn();
    $('#gd-analyze-variant-filter-container').fadeIn();

    // Recreate the initial DataTable.
    this.drawDatatable();

    // Kick off request to update the Variant data being displayed.
    this.updateVariantList();
  },


  /** Redraws the datatable based on the selection. */
  drawDatatable: function() {
    this.datatableComponent = new gd.ServerSideDataTableComponent({
        el: $('#gd-datatable-hook'),
        serverTarget: '/_/variants',
        serverParamsInjector: _.bind(
            this.addFilterDataToDataTableServerSideRequest, this)
    });

    // Listen to events from the DataTable to provide appropriate UI affordance.
    this.listenTo(this.datatableComponent, 'START_LOADING',
        _.bind(this.setUIStartLoadingState, this));
    this.listenTo(this.datatableComponent, 'DONE_LOADING',
        _.bind(this.setUIDoneLoadingState, this));
  },


  /** Kicks off the process for updating the list of displayed variants. */
  updateVariantList: function() {
    // Update the UI to show that the new Variant list is loading.
    this.setUIStartLoadingState();

    var requestData = {
      'projectUid': this.model.get('project').uid,
      'refGenomeUid': this.model.get('refGenomeUid'),
      'variantFilterString': $('#gd-new-filter-input').val(),
      'melt': $('input:radio[name=melt]:checked').val()
    };

    if (this.visibleKeyNames) {
      requestData['visibleKeyNames'] = JSON.stringify(this.visibleKeyNames);
    }

    $.get('/_/variants', requestData,
        _.bind(this.handleGetVariantListResponse, this));
  },


  /** Handles the server response containing Variants data. */
  handleGetVariantListResponse: function(response) {
    // Parse Variant data.
    var variantListData = JSON.parse(response.variant_list_json);
    if (variantListData.obj_list.length) {
      this.variantList = variantListData.obj_list;
      this.fieldConfig = variantListData.field_config;
    } else {
      this.variantList = [];
      this.fieldConfig = {};
    }

    var numTotalVariants = Number(response.num_total_variants);

    // Parse VariantSet data.
    this.variantSetList = JSON.parse(response.variant_set_list_json).obj_list;

    // Grab the field key map data.
    this.variantKeyMap = JSON.parse(response.variant_key_filter_map_json);

    // Redraw the datatable.
    this.datatableComponent.update(this.variantList, this.fieldConfig,
        numTotalVariants);

    // Does this need to be called every time?
    this.drawDropdowns();

    // Reset the ui.
    this.setUIDoneLoadingState();
  },


  /** Handles a click on the 'Apply Filter' button. */
  handleApplyFilterClick: function() {
    // Kick off request to update the Variant data being displayed.
    this.updateVariantList();
  },


  /** Provide affordance while doing Variant query. */
  setUIStartLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 0.5);
  },


  /** Reset UI changes after loading complete.. */
  setUIDoneLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 1);
  },


  /**
   * Updates the array of data representing GET params to be sent to the
   * server by DataTables. See jquery.dataTables.js.
   *
   * @param {Array.<Object.<string, string>>} aoData Array of objects with keys
   *     'name' and 'value'. The array is passed by reference so this method
   *      must mutate this array.
   */
  addFilterDataToDataTableServerSideRequest: function(aoData) {
    var requestData = {
      'projectUid': this.model.get('project').uid,
      'refGenomeUid': this.model.get('refGenomeUid'),
      'variantFilterString': $('#gd-new-filter-input').val(),
      'melt': $('input:radio[name=melt]:checked').val()
    };

    _.each(_.pairs(requestData), function(pair) {
      aoData.push({'name': pair[0], 'value': pair[1]});
    })
  },


  /* Adds options to the variant set dropdown menu. */
  drawDropdowns: function() {
    // Create 'add selected to set' button, display list of available sets
    this.addVariantSetDropdownSubmenu(
        this.variantSetList, 'add', 'Add to set');
    this.addVariantSetDropdownSubmenu(
        this.variantSetList, 'remove', 'Remove from set');
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

    this.datatableComponent.addDropdownOption(setListHTML, '');
  },


  /**
   * Create a single html list object for a sub-dropdown set option.
   *
   * @param variantSet The variant set dictionary object from the
   *     django adapter.
   * @param variantSetAction The action to perform, either 'add' or 'remove'
   */
  addIndividualSetListItem: function(variantSet, variantSetAction) {
    return (
      '<li><a tabindex="-1" class="gd-variant-set-action" ' +
        'data-variant-set-uid="' + variantSet.uid + '" ' +
        'data-variant-set-action="' + variantSetAction +
        '" href="#">' + variantSet.label + '</a></li>');
  },


  /** Handles a click on a variant set action button. */
  handleVariantSetClick: function(ev) {
    // Post to this same view for now.
    var postUrl = '/_/variants/modify_set_membership';

    // Grab the selected rows.
    var postData = {
        projectUid: this.model.get('project').uid,
        refGenomeUid: this.model.get('refGenomeUid'),
        variantUidList: this.datatableComponent.getCheckedRowUids(),
        variantSetAction: $(ev.target).data('variant-set-action'),
        variantSetUid: $(ev.target).data('variant-set-uid')
    };

    var validationResult = this.validatePostData(postData);
    if (!validationResult.is_success) {
      this.showAlertMessage(validationResult.error_msg, 'error');
      return;
    }

    var onSuccess = _.bind(function(response) {
      // Repopulate json data from response
      this.variantList = JSON.parse(response.variant_list_json).obj_list;
      this.variantSetList = JSON.parse(response.variant_set_list_json).obj_list;

      // Show alert message based on response.
      this.showAlertMessage(response.alert_msg, response.alert_type);

      // Redraw the datatable.
      this.datatableComponent.update(this.variantList);
      this.drawDropdowns();
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
  },


  /** Shows an alert message. */
  showAlertMessage: function(msg, alertType) {
    $('#gd-variant-set-action-submit-alert-msg').text(msg);
    $('#gd-variant-set-action-submit-alert').removeClass(
        'alert-warn alert-info alert-error');
    $('#gd-variant-set-action-submit-alert').addClass('alert-' + alertType);
    $('#gd-variant-set-action-submit-alert').show();
  },


  /**
   * Creates the modal for selecting which fields are being displayed and
   * shows it.
   */
  handleShowFieldSelect: function() {
    var variantFieldSelectModel = new Backbone.Model({
        'variantKeyMap': this.variantKeyMap
    });
    this.variantFieldSelect = new gd.VariantFieldSelectComponent(
        {'model': variantFieldSelectModel});
    this.listenTo(this.variantFieldSelect, 'updateSelectedFields',
        _.bind(this.handleUpdateSelectedFields, this));
  },


  /** Handles an update event from the field select component. */
  handleUpdateSelectedFields: function(selectedKeyNames) {
    // Hide the modal.
    this.variantFieldSelect.hide();

    // Update the attribute that stores the visible field names.
    this.visibleKeyNames = selectedKeyNames;

    // Kick off a request to update the view.
    this.updateVariantList();
  }
},

/** Static properties */
{
  /**
   * Default number of items to show in the DataTable at one time.
   * @type {number}
   */
  DEFAULT_PAGE_SIZE: 100
}
);
