/**
 * Component for modal that allows uploading a single sample through the
 * browser.
 */


gd.UploadSingleSampleModal= Backbone.View.extend({

  initialize: function() {
    if (!this.model) {
      throw "gd.UploadSingleSampleModal requires model.";
      return;
    }

    this.render();
  },

  render: function() {
    this.listenToControls();
  },

  listenToControls: function() {
    $('#gd-samples-upload-single-submit').click(
        _.bind(this.handleUploadThroughBrowser, this));
  },

  enterLoadingState: function() {
    $(".gd-samples-upload-single-submit")
        .prop('disabled', true);

    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin();
  },

  exitLoadingState: function() {
    $(".gd-samples-upload-single-submit")
        .prop('disabled', false);
    gd.Spinner.globalClear();
  },

  /** Uploads genome through browser. */
  handleUploadThroughBrowser: function() {
    // Parse the inputs.
    var requestData = this.prepareRequestData(
        'gd-samples-upload-single-form');

    // Validate the request client-side, even though we just grab
    // the data again into FormData form below.
    if (!this.validateUploadThroughBrowserForm(requestData)) {
      return;
    }

    // Put UI in loading state.
    this.enterLoadingState();

    var onSuccess = _.bind(function(responseData) {
      this.exitLoadingState();

      $('.modal').modal('hide');
    }, this);

    var formData = new FormData(
        $('#gd-samples-upload-single-form')[0]);
    $.ajax({
      url: '/_/samples/upload_single_sample',
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

  /** Parses the form files and prepares the data. */
  prepareRequestData: function(formId) {
    var requestData = {}

    // This works correctly for all inputs except radio button.
    // See below for that.
    var formInputs = $('#' + formId + ' :input');
    _.each(formInputs, function(inputObj) {
      requestData[inputObj.name] = inputObj.value;
    });

    return requestData;
  },


  validateUploadThroughBrowserForm: function(requestData) {
    if (!requestData['sampleLabel'].length) {
      alert('Please enter a name for this sample.');
      return false;
    }

    if (!requestData['fastq1'].length) {
      alert('Please specify a file location.');
      return false;
    }

    return true;
  },

  destroy: function() {
    // TODO(gleb): Is there a simpler way to unbind all events?
    $('#gd-samples-upload-single-submit').unbind('click');
    this.remove();
    this.unbind();
  }
});
