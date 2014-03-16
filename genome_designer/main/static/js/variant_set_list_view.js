/**
 * @fileoverview List of variant sets.
 */


gd.VariantSetListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-variant-sets').addClass('active');

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-variant_set_list_view-datatable-hook'),
        serverTarget: '/_/sets',
        controlsTemplate: '/_/templates/variant_set_list_controls',
        requestData: {projectUid: this.model.get('uid')}
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.listenToControls, this));
  },

  listenToControls: function() {
    $('#gd-variant-set-form-from-file-submit').click(
        _.bind(this.handleFormSubmitFromFile, this));
    $('#gd-variant-set-form-empty-submit').click(
        _.bind(this.handleFormSubmitEmpty, this));
  },

  /**
   * Handles submitting the request to create a new VariantSet from a file
   * containing a list of Variants.
   */
  handleFormSubmitFromFile: function() {
    // We validate by parsing the input elements, but then use the HTML5
    // FormData API to actually send the data to the server. FormData
    // makes the file upload "just work".
    var formDataForValidation = this.prepareRequestData(
        'gd-variant-set-create-form-from-file');
    if (!this.validateCreateFromFileRequestData(formDataForValidation)) {
      return;
    }

    var onSuccess = function(responseData) {
      if (responseData.error.length) {
        alert('Error creating variant set: ' + responseData.error);
        return;
      }

      // Success, reload the page.
      window.location.reload();
    }

    var formData = new FormData($('#gd-variant-set-create-form-from-file')[0]);
    $.ajax({
      url: '/_/sets/create',
      type: 'POST',
      data: formData,
      success: onSuccess,

      // The following 3 param settings are necessary for properly passing
      // formData. See: http://stackoverflow.com/questions/166221/how-can-i-upload-files-asynchronously-with-jquery
      cache: false,
      contentType: false,
      processData: false
    });
  },

  /** Handles creating a new empty Variant set. */
  handleFormSubmitEmpty: function() {
    // Parse the inputs.
    var requestData = this.prepareRequestData('gd-variant-set-form-empty');

    // Validate the request client-side.
    if (!this.validateCreateEmptyRequestData(requestData)) {
      return;
    }

    // Make the request.
    $.post('/_/sets/create', requestData, function(responseData) {
      // Check for error and show in ui. Don't reload the page.
      if (responseData.error.length) {
        alert('Error creating variant set: ' + responseData.error);
        return;
      }

      // Success, reload the page.
      window.location.reload();
    });
  },

  /** Parses the form files and prepares the data. */
  prepareRequestData: function(formId) {
    var requestData = {}
    var formInputs = $('#' + formId + ' :input');
    _.each(formInputs, function(inputObj) {
      requestData[inputObj.name] = inputObj.value;
    });
    return requestData;
  },

  /** Validates the request to be sent to the server. */
  validateCreateEmptyRequestData: function(requestData) {
    if (!this.validateCommon(requestData)) {
      return false;
    }

    return true;
  },

  /** Validates the request to be sent to the server. */
  validateCreateFromFileRequestData: function(requestData) {
    if (!this.validateCommon(requestData)) {
      return false;
    }

    if (!requestData['vcfFile'].length) {
      alert('Please select a vcf file to upload.');
      return false;
    }

    return true;
  },

  /** Validation common to both forms. */
  validateCommon: function(requestData) {
    if (!requestData['refGenomeUid'].length) {
      alert('Please select a reference genome.');
      return false;
    }

    if (!requestData['variantSetName'].length) {
      alert('Please enter a variant set name.');
      return false;
    }

    return true;
  }
});
