/** * @fileoverview Sub-view displaying Variants.
 */


gd.TabAnalyzeSubviewVariants = gd.TabAnalyzeSubviewAbstractBase.extend(
/** Prototype properties */
{
  /** @override */
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

    // Check whether a filter is set on the view.
    var paramObject = gd.Util.getQueryStringAsObject();
    var maybeQueryString = window.location.search;
    if ('filter' in paramObject) {
      var filterString = paramObject['filter'];
      this.model.set('filterString', filterString);
    }

    gd.TabAnalyzeSubviewAbstractBase.prototype.initialize.call(this);
  },


  /** @override */
  renderControls: function() {
    // Request the filter control html from the server for now.
    // TODO: Client-side template support.
    $.get('/_/templates/variant_filter_controls', _.bind(function(response) {
      // Append the DOM.
      $('#gd-analyze-subview-controls-hook').append(response);

      // Fill in the filter if relevant.
      $('#gd-new-filter-input').val(this.model.get('filterString'));

      // Register listeners.
      $('#gd-filter-box-apply-btn').click(
          _.bind(this.handleApplyFilterClick, this));
      $('#gd-filter-field-select-btn').click(
          _.bind(this.handleShowFieldSelect, this));
      $('#gd-filter-refresh-materialized').click(
          _.bind(this.handleRefreshMaterializedView, this));
    }, this));
  },


  /** @override */
  renderDatatable: function() {
    // Recreate the initial DataTable.
    this.createDatatableComponent();

    // Kick off request to update the Variant data being displayed.
    this.updateVariantList();
  },


  /** Redraws the datatable based on the selection. */
  createDatatableComponent: function() {
    this.datatableComponent = new gd.ServerSideDataTableComponent({
        el: $('#gd-datatable-hook'),
        serverTarget: '/_/variants',
        serverParamsInjector: _.bind(
            this.addFilterDataServerSideRequest, this)
    });

    // Listen to events from the DataTable to provide appropriate UI affordance.
    this.listenTo(this.datatableComponent, 'START_LOADING',
        _.bind(this.setUIStartLoadingState, this));
    this.listenTo(this.datatableComponent, 'DONE_LOADING',
        _.bind(this.setUIDoneLoadingState, this));
  },


  /**
   * Updates the array of data representing GET params to be sent to the
   * server by DataTables. See jquery.dataTables.js.
   *
   * @param {Array.<Object.<string, string>>} aoData Array of objects with keys
   *     'name' and 'value'. The array is passed by reference so this method
   *      must mutate this array.
   */
  addFilterDataServerSideRequest: function(aoData) {
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


  /** Provide affordance while doing Variant query. */
  setUIStartLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 0.5);
  },


  /** Reset UI changes after loading complete.. */
  setUIDoneLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 1);
  },


  /** Handles the server response containing Variants data. */
  handleGetVariantListResponse: function(response) {
    if ('error' in response) {
      this.handleGetVariantListError(response.error);
      return;
    }

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


  handleGetVariantListError: function(errorMsg) {
    $('#gd-snp-filter-error-msg').text(errorMsg);
    $('#gd-snp-filter-error').show();
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


  /**
   * Sends a request to the server to refresh the materialized data table and
   * reloads the page on success.
   */
  handleRefreshMaterializedView: function() {
    this.setUIStartLoadingState();

    var requestData = {
      'refGenomeUid': this.model.get('refGenomeUid'),
    }
    $.get('/_/variants/refresh_materialized_variant_table', requestData,
        function(data) {
          window.location.reload();
        });
  },


  /** Handles an update event from the field select component. */
  handleUpdateSelectedFields: function(selectedKeyNames) {
    // Hide the modal.
    this.variantFieldSelect.hide();

    // Update the attribute that stores the visible field names.
    this.visibleKeyNames = selectedKeyNames;

    // Kick off a request to update the view.
    this.updateVariantList();
  },


  /** Kicks off the process for updating the list of displayed variants. */
  updateVariantList: function() {
    // Update the UI to show that the new Variant list is loading.
    this.setUIStartLoadingState();

    var requestData = {
      'projectUid': this.model.get('project').uid,
      'refGenomeUid': this.model.get('refGenomeUid'),
      'variantFilterString': this.model.get('filterString'),
      'melt': $('input:radio[name=melt]:checked').val()
    };

    if (this.visibleKeyNames) {
      requestData['visibleKeyNames'] = JSON.stringify(this.visibleKeyNames);
    }

    $.get('/_/variants', requestData,
        _.bind(this.handleGetVariantListResponse, this));
  },


  /** Handles a click on the 'Apply Filter' button. */
  handleApplyFilterClick: function() {
    // Update the model.
    // NOTE: We are not being super consistent about this but can't hurt
    //     to keep doing it where possible.
    this.model.set('filterString', $('#gd-new-filter-input').val());

    // Update the url with the filter.
    var updatedPath = gd.Util.updateQueryStringParameter(
        gd.Util.getFullPath(), 'filter', this.model.get('filterString'));
    this.trigger('NAVIGATE', {'path': updatedPath});

    // Kick off request to update the Variant data being displayed.
    this.updateVariantList();
  },


  /* Adds options to the variant set dropdown menu. */
  drawDropdowns: function() {
    // Create 'add selected to set' button, display list of available sets
    this.addVariantSetDropdownSubmenu(
        this.variantSetList, 'add', 'Add to set');
    this.addVariantSetDropdownSubmenu(
        this.variantSetList, 'remove', 'Remove from set');

    $('.gd-id-variant-set-action').click(
        _.bind(this.handleVariantSetActionClick, this));
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
      '<li><a class="gd-id-variant-set-action" ' +
          'data-variant-set-uid="' + variantSet.uid + '" ' +
          'data-variant-set-action="' + variantSetAction +
          '" href="#">' +
        variantSet.label +
      '</a></li>');
  },


  /** Handles a click on a variant set action button. */
  handleVariantSetActionClick: function(ev) {
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

    var onSuccess = _.bind(function(response_json) {
      var response = JSON.parse(response_json);
      if (response.alert_type == 'error') {
        alert(response.alert_msg);
        return;
      }

      this.handleApplyFilterClick();
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


  /**
   * Shows an alert message.
   *
   * TODO: Add back the DOM hooks to make this work.
   */
  showAlertMessage: function(msg, alertType) {
    $('#gd-variant-set-action-submit-alert-msg').text(msg);
    $('#gd-variant-set-action-submit-alert').removeClass(
        'alert-warn alert-info alert-error');
    $('#gd-variant-set-action-submit-alert').addClass('alert-' + alertType);
    $('#gd-variant-set-action-submit-alert').show();
  },


  /** Properly clean up this component. Inheriting classes should override. */
  destroy: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy()
    }
    this.datatableComponent = null;
  }
},

/** Static properties */
{
  /**
   * Default number of items to show in the DataTable at one time.
   * @type {number}
   */
  DEFAULT_PAGE_SIZE: 100,
}
);
