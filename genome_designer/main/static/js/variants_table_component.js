/**
 * @fileoverview Component that shows Variants.
 */


gd.VariantsTableComponent = Backbone.View.extend({

  /** Override. */
  initialize: function() {
    // Populate model with attributes.
    this.model.set('filterString', 'filterString' in this.options ?
        this.options['filterString'] : '');

    this.model.set('filterDisabled', 'filterDisabled' in this.options ?
        this.options['filterDisabled'] : false);

    this.model.set('isMelted', false);

    this.render();
  },

  /** Override. */
  render: function() {
    this.renderDatatable();

  },

  /** @override */
  renderDatatable: function() {
    // Always start this flow by checking whether the materialized view
    // is valid. Then do the rest of the create after.
    this.prepareMaterializedView(
        _.bind(this.afterMaterializedViewReady, this));
  },

  /** @override */
  decorateControls: function() {
    var filterEl = $('#gd-new-filter-input');
    var applyFilterBtnEl = $('#gd-filter-box-apply-btn');

    // Decorate the filter element.
    filterEl.val(this.model.get('filterString'));
    if (this.model.get('filterDisabled')) {
      filterEl.prop('disabled', true);
      applyFilterBtnEl.hide();
    }

    // Decorate the cast/melt toggle.
    if (this.model.get('isMelted')) {
      $('#gd-snp-filter-melt-toggle-melt').addClass('active');
    } else {
      $('#gd-snp-filter-melt-toggle-cast').addClass('active');
    }
    $('.gd-snp-filter-melt-toggle').click(
        _.bind(this.handleMeltedToggleClick, this));

    // Manually handle error alert close.
    $('#gd-snp-filter-error-close-btn').click(
        _.bind(this.handleErrorAlertClose, this));
  },

  afterMaterializedViewReady: function() {
    // Recreate the initial DataTable.
    this.createDatatableComponent();

    // // Kick off request to update the Variant data being displayed.
    this.updateVariantList();
  },

  /**
   * Checks if the underlying materialized view needs to be updated
   * and kick off a whole ui affordance to manage that.
   */
  prepareMaterializedView: function(callback) {
    // First make a request to check if needs to be refreshed.
    var requestData = {
      'refGenomeUid': this.model.get('refGenomeUid'),
    };
    $.get('/_/variants/is_materialized_view_valid', requestData,
        _.bind(function(responseData) {
          if (responseData.isValid) {
            callback.call();
          } else {
            this.handleRefreshMaterializedView(callback);
          }
        }, this));
  },

  /** Redraws the datatable based on the selection. */
  createDatatableComponent: function() {
    this.datatableComponent = new gd.ServerSideDataTableComponent({
        el: $('#gd-datatable-hook'),
        serverTarget: '/_/variants',
        controlsTemplate: '/_/templates/variant_filter_controls',
        requestData: {'refGenomeUid': this.model.get('refGenomeUid')}
    });

    // Listen to events from the DataTable to provide appropriate UI affordance.
    this.listenTo(this.datatableComponent, 'START_LOADING',
        _.bind(this.setUIStartLoadingState, this));
    this.listenTo(this.datatableComponent, 'DONE_LOADING',
        _.bind(this.setUIDoneLoadingState, this));
    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  /** Provide affordance while doing Variant query. */
  setUIStartLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 0.5);

    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin()
  },

  /** Reset UI changes after loading complete.. */
  setUIDoneLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 1);
    if (this.loadingSpinner) {
      this.loadingSpinner.stop();
    }
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

    this.datatableComponent.update(this.variantList, this.fieldConfig,
        numTotalVariants);

    // Reset the ui.
    this.setUIDoneLoadingState();
  },

  handleGetVariantListError: function(errorMsg) {
    this.setUIDoneLoadingState();
    $('#gd-snp-filter-error-msg').text(errorMsg);
    $('#gd-snp-filter-error').show();
  },

  /**
   * Hides the alert box.
   *
   * Previously, we were adding Bootstrap's data-dismiss="alert" to the dom
   * but that was causing the entire element to be deleted so that subsequent
   * alerts wouldn't show.
   */
  handleErrorAlertClose: function() {
    $('#gd-snp-filter-error').hide();
  },

  /**
   * Sends a request to the server to refresh the materialized data table and
   * reloads the page on success.
   */
  handleRefreshMaterializedView: function(onSuccess) {
    this.setUIStartLoadingState();

    var requestData = {
      'refGenomeUid': this.model.get('refGenomeUid'),
    };
    $.get('/_/variants/refresh_materialized_variant_table', requestData,
        _.bind(function(data) {
          this.setUIDoneLoadingState();
          onSuccess();
        }, this));
  },

  /** Kicks off the process for updating the list of displayed variants. */
  updateVariantList: function() {
    // Update the url with the filter.
    var updatedPath = gd.Util.updateQueryStringParameter(
        gd.Util.getFullPath(), 'filter', this.model.get('filterString'));
    updatedPath = gd.Util.updateQueryStringParameter(
        updatedPath, 'melt', this.model.get('is_melted') ? '1' : '0');
    this.trigger('NAVIGATE', {'path': updatedPath});

    // Hide any alert.
    this.handleErrorAlertClose();

    // Update the UI to show that the new Variant list is loading.
    this.setUIStartLoadingState();

    var requestData = {
      'projectUid': this.model.get('project').uid,
      'refGenomeUid': this.model.get('refGenomeUid'),
      'variantFilterString': this.model.get('filterString'),
      'melt': this.model.get('isMelted') ? 1 : 0
    };

    if (this.visibleKeyNames) {
      requestData['visibleKeyNames'] = JSON.stringify(this.visibleKeyNames);
    }

    $.get('/_/variants', requestData,
        _.bind(this.handleGetVariantListResponse, this));
  },

  /** Handles a click on one of the melted toggle buttons. */
  handleMeltedToggleClick: function(e) {
    // Update the model.
    this.model.set('isMelted', Boolean($(e.target).data('melted')));
    this.model.set('filterString', $('#gd-new-filter-input').val());
    this.updateVariantList();
  },

  /** Clean up the component. */
  destroy: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy()
    }
    this.datatableComponent = null;
  }
});
