/**
 * @fileoverview Sub-view displaying Contigs.
 */


gd.TabAnalyzeSubviewContigs = gd.TabAnalyzeSubviewAbstractBase.extend({

  /** @override */
  initialize: function() {
    this.render();
  },

  render: function() {
    this.createDatatableComponent();

    // Kick off request to update the Contig data being displayed.
    this.updateContigList();
  },

  /** Draws or redraws the table. */
  createDatatableComponent: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    var requestData = {
        alignmentGroupUid: this.model.get('alignmentGroupUid')
    };

    this.datatableComponent = new gd.ServerSideDataTableComponent({
        el: $('#gd-datatable-hook-contigs'),
        serverTarget: '/_/contigs',
        controlsTemplate: '/_/templates/contig_list_controls',
        requestData: requestData,
        serverParamsInjector: _.bind(
            this.addFilterDataServerSideRequest, this)
    });

    this.listenToOnce(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  /** Returns request data object for making a request for variants. */
  _prepareContigListRequestData: function() {
    var requestData = {
      'alignmentGroupUid': this.model.get('alignmentGroupUid')
    };
    return requestData;
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
    var requestData = this._prepareContigListRequestData();

    _.each(_.pairs(requestData), function(pair) {
      aoData.push({'name': pair[0], 'value': pair[1]});
    })
  },

  /** Kicks off the process for updating the list of displayed contigs. */
  updateContigList: function() {
    // Update the UI to show that the new Contig list is loading.
    this.setUIStartLoadingState();

    var requestData = this._prepareContigListRequestData();

    $.get('/_/contigs', requestData,
        _.bind(this.handleGetContigListResponse, this));
  },

  /** Handles the server response containing Contig data. */
  handleGetContigListResponse: function(response) {
    // Parse Contig data.
    var contigListData = JSON.parse(response.contig_list_json);
    if (contigListData.obj_list.length) {
      this.contigList = contigListData.obj_list;
      this.fieldConfig = contigListData.field_config;
    } else {
      this.contigList = [];
      this.fieldConfig = {};
    }

    var numTotalVariants = Number(response.num_total_variants);

    var timeForLastResult = .68;

    this.datatableComponent.update(this.contigList, this.fieldConfig,
        numTotalVariants, timeForLastResult);

    // Reset the ui.
    this.setUIDoneLoadingState();
  },

  /** Provide affordance while doing Variant query. */
  setUIStartLoadingState: function() {
    $('#gd-datatable-hook-datatable_wrapper').css('opacity', 0.5);

    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin()
  },

  /** Reset UI changes after loading complete.. */
  setUIDoneLoadingState: function() {
    $('#gd-datatable-hook-contigs-datatable_wrapper').css('opacity', 1);
    gd.Spinner.globalClear();
  },

  decorateControls: function() {
    if (this.contigControlsComponent) {
      this.contigControlsComponent.destroy();
    }

    this.contigControlsComponent = new gd.ContigControlsComponent({
      el: '#gd-datatable-hook-control',
      datatableComponent: this.datatableComponent,
      alignmentGroupUid: this.model.attributes.alignmentGroupUid
    });

    this.listenToOnce(this.contigControlsComponent, 'MODELS_UPDATED',
        _.bind(this.render, this));
  }
});
