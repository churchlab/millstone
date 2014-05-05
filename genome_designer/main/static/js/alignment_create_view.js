/**
 * @fileoverview List of variants.
 */


gd.AlignmentCreateView = Backbone.View.extend({
  el: '#gd-page-container',

  /** Override. */
  initialize: function() {
    this.render();
  },

  /** Override. */
  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.nameInput = $('#gd-alignment-create-name-input');

    this.redrawRefGenomeDatatable();
    this.redrawSampleControlsDatatable();
  },

  /** Draws or redraws the table. */
  redrawRefGenomeDatatable: function() {
    if (this.refGenomeDataTable) {
      this.refGenomeDataTable.destroy();
    }

    this.refGenomeDataTable = new gd.DataTableComponent({
        el: $('#gd-datatable-ref_genome-hook'),
        serverTarget: '/_/ref_genomes',
        controlsTemplate: '/_/templates/reference_genome_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });
    this.listenTo(this.refGenomeDataTable, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateRefGenomeControls, this));
  },

  /** Create component to wrap ref genome controls, and listen for events. */
  decorateRefGenomeControls: function() {
    this.refGenomeControlsComponent = new gd.RefGenomeControlsComponent({
      el: '#gd-datatable-ref_genome-hook-control'
    });

    this.listenTo(this.refGenomeControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawRefGenomeDatatable, this));
  },

  /** Draws or redraws the table. */
  redrawSampleControlsDatatable: function() {
    if (this.samplesDatatable) {
      this.samplesDatatable.destroy();
    }

    this.samplesDatatable = new gd.DataTableComponent({
        el: $('#gd-datatable-samples-hook'),
        serverTarget: '/_/samples',
        controlsTemplate: '/_/templates/sample_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });

    this.listenTo(this.samplesDatatable, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateSamplesControls, this));
  },

  /** Create component to wrap samples controls, and listen for events. */
  decorateSamplesControls: function() {
    this.samplesControlsComponent = new gd.SamplesControlsComponent({
      el: '#gd-datatable-samples-hook-control',
      model: this.model,
      datatableComponent: this.samplesDatatable
    });

    this.listenTo(this.samplesControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawSampleControlsDatatable, this));
  },

  events: {
    'click #gd-align-create-submit-btn': 'handleSubmitClick',
    'click #gd-align-create-submit-error-close': 'handleErrorAlertClose'
  },

  /**
   * Handles a click on the submit button, aggregating the config
   * options selected and sending a request to the server to start
   * an alignment run.
   */
  handleSubmitClick: function() {
    // Disable the button to let the user know things are happening.
    $('#gd-align-create-submit-btn').removeClass("btn-success");
    $('#gd-align-create-submit-btn').addClass("disabled");

    // TODO: Show a spinner.

    // Post to this same view for now.
    var postUrl = window.location.pathname;

    // Grab the selected rows.
    var postData = {
        name: this.nameInput.val(),
        refGenomeUidList: this.refGenomeDataTable.getCheckedRowUids(),
        sampleUidList: this.samplesDatatable.getCheckedRowUids()
    };

    // Validate the data.
    var validationResult = this.validatePostData(postData);
    if (!validationResult.is_success) {
      $('#gd-align-create-submit-error-msg').text(validationResult.error_msg);
      $('#gd-align-create-submit-error').show();
      $('#gd-align-create-submit-btn').addClass("btn-success");
      $('#gd-align-create-submit-btn').removeClass("disabled");
      return;
    }

    var onSuccess = function(data) {
      if ('error' in data) {
        $('#gd-align-create-submit-error-msg').text(data.error);
        $('#gd-align-create-submit-error').show();
        $('#gd-align-create-submit-btn').addClass("btn-success");
        $('#gd-align-create-submit-btn').removeClass("disabled");
        return;
      }

      // Redirect according to the server response.
      window.location.href = data.redirect;
    };

    // Execute the post. Handle the response accordingly.
    $.post(postUrl, JSON.stringify(postData), onSuccess, 'json');
  },

  /**
   * Validate the post data before submitting.
   * @param {object} postData The creation data to be posted to the server.
   * @return {object} Response with keys:
   *     * is_success {boolean} Whether validation succeded.
   *     * error_msg {string} Human-readable description of the error.
   */
  validatePostData: function(postData) {
    if (!postData.name) {
      return {
          is_success: false,
          error_msg: 'Name required.'
      };
    }

    if (!postData.refGenomeUidList.length) {
      return {
          is_success: false,
          error_msg: 'Please select at least one reference genome.'
      };
    }

    if (!postData.sampleUidList.length) {
      return {
          is_success: false,
          error_msg: 'Please select at least one sample.'
      };
    }

    return {
          is_success: true,
          error_msg: ''
    };
  },

  /**
   * Hides the alert box.
   *
   * Previously, we were adding Bootstrap's data-dismiss="alert" to the dom
   * but that was causing the entire element to be deleted so that subsequent
   * alerts wouldn't show.
   */
  handleErrorAlertClose: function() {
    $('#gd-align-create-submit-error').hide();
  },
});
