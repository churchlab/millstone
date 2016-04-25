/**
 * @fileoverview Component that decorates the controls for the list of
 *     Variant Sets.
 */


gd.VariantsetsControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);

    if (!this.model) {
      throw "VariantSetControlsComponent requires model.";
    }
    this.render();
  },

  render: function() {
    this.decorateControls();
  },

  decorateControls: function() {
    this.drawDropdownOptions();
    this.listenToControls();
  },

  listenToControls: function() {
    $('#gd-variant-set-form-from-file-submit').click(
        _.bind(this.handleFormSubmitFromFile, this));
    $('#gd-variant-set-form-empty-submit').click(
        _.bind(this.handleFormSubmitEmpty, this));
  },

  /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete samples.
    var deleteOptionHtml =
        '<a href="#" class="gd-id-variant-sets-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-variant-sets-delete').click(_.bind(this.handleDelete, this));
  },

  /** Sends request to delete selected samples. */
  handleDelete: function() {
    var variantSetUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!variantSetUidList.length) {
      alert("Please select variant sets to delete.");
      return;
    }

    // Get confirmation from user.
    var agree = confirm("Are you sure you want delete these variant sets?");
    if (!agree) {
      return;
    }

    this.setUIStartLoadingState();

    var postData = {
        variantSetUidList: variantSetUidList,
    };

    $.post('/_/variants/delete', JSON.stringify(postData),
        _.bind(this.handleDeleteResponse, this));
  },

  handleDeleteResponse: function(response) {
    this.setUIDoneLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
  },

  /** Show loading feedback while loading. */
  setUIStartLoadingState: function() {
    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin()
  },

  /** Reset UI changes after loading complete.. */
  setUIDoneLoadingState: function() {
    gd.Spinner.globalClear();
  },

  /**
   * Handles submitting the request to create a new VariantSet from a file
   * containing a list of Variants.
   */
  handleFormSubmitFromFile: function() {
    var submitButtonEl = $('#gd-variant-set-form-from-file-submit');
    submitButtonEl.prop('disabled', true);

    // We validate by parsing the input elements, but then use the HTML5
    // FormData API to actually send the data to the server. FormData
    // makes the file upload "just work".
    var formDataForValidation = this.prepareRequestData(
        'gd-variant-set-create-form-from-file');
    if (!this.validateCreateFromFileRequestData(formDataForValidation)) {
      submitButtonEl.prop('disabled', false);
      return;
    }

    var onSuccess = function(responseData) {
      if ('error' in responseData && responseData.error.length) {
        alert('Error creating variant set: ' + responseData.error);
        submitButtonEl.prop('disabled', false);
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
    var submitButtonEl = $('#gd-variant-set-form-empty-submit');
    submitButtonEl.prop('disabled', true);

    // Parse the inputs.
    var requestData = this.prepareRequestData('gd-variant-set-form-empty');

    // Validate the request client-side.
    if (!this.validateCreateEmptyRequestData(requestData)) {
      submitButtonEl.prop('disabled', false);
      return;
    }

    // Make the request.
    $.post('/_/sets/create', requestData, function(responseData) {
      // Check for error and show in ui. Don't reload the page.
      if (responseData.error.length) {
        alert('Error creating variant set: ' + responseData.error);
        submitButtonEl.prop('disabled', false);
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
  },

  destroy: function() {
    $('#gd-variant-set-form-from-file-submit').unbind('click')
    $('#gd-variant-set-form-empty-submit').unbind('click')
    this.remove();
    this.unbind();
  }
});
