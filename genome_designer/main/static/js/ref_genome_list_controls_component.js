/**
 * @fileoverview Component that decorates the controls for the list of
 *     ReferenceGenomes.
 */


gd.RefGenomeControlsComponent = Backbone.View.extend({
  initialize: function() {
    this.listenToControls();

    this.maybeCreateFineS3Uploader();
  },

  listenToControls: function() {
    $('#gd-ref-genome-upload-from-server-location-submit').click(
        _.bind(this.handleCreateFromServerLocation, this));
    $('#gd-ref-genome-create-from-ncbi-submit').click(
        _.bind(this.handleCreateFromNCBI, this));
  },

  /** Puts UI in the loading state. */
  enterLoadingState: function() {
    $(".gd-id-form-submit-button")
        .prop('disabled', true);

    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin();
  },

  /** Puts UI in the loading state. */
  exitLoadingState: function() {
    $(".gd-id-form-submit-button")
        .prop('disabled', false);
    this.loadingSpinner.stop();
  },

  handleCreateFromServerLocation: function() {
    // Parse the inputs.
    var requestData = this.prepareRequestData(
        'gd-ref-genome-upload-from-server-location');

    // Validate the request client-side.
    if (!this.validateCreateFromServerLocation(requestData)) {
      return;
    }

    // Put UI in loading state.
    this.enterLoadingState();

    // Make the request.
    $.post('/_/ref_genomes/create_from_server_location', requestData,
        _.bind(function(responseData) {
          // Check for error and show in ui. Don't reload the page.
          if (responseData.error.length) {
            alert('Error creating reference genome: ' + responseData.error);
            this.exitLoadingState();
            return;
          }

          // Success, reload the page.
          window.location.reload();
        }, this));
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

    // The above screws up the radio button. Do it manually here.
    requestData['importFileFormat'] =
        $('input[name="importFileFormat"]:checked').val();

    return requestData;
  },

  /** Validation common to both forms. */
  validateCreateFromServerLocation: function(requestData) {
    if (!requestData['refGenomeLabel'].length) {
      alert('Please enter a name for this reference genome.');
      return false;
    }

    if (!requestData['refGenomeFileLocation'].length) {
      alert('Please specify a file location.');
      return false;
    }

  if (!requestData['importFileFormat'] in ['genbank', 'fasta']) {
      alert('Please select format for file.');
      return false;
    }

    return true;
  },

  handleCreateFromNCBI: function() {
    // Parse the inputs.
    var requestData = this.prepareRequestData(
        'gd-ref-genome-create-from-ncbi');

    // Validate the request client-side.
    if (!this.validateCreateFromNCBILocation(requestData)) {
      return;
    }

    // Put UI in loading state.
    this.enterLoadingState();

    // Make the request.
    $.post('/_/ref_genomes/create_from_ncbi', requestData,
        _.bind(function(responseData) {
          // Check for error and show in ui. Don't reload the page.
          if (responseData.error.length) {
            alert('Error creating reference genome: ' + responseData.error);
            this.exitLoadingState();
            return;
          }

          // Success, reload the page.
          window.location.reload();
        }, this));
  },

  /** Validation common to both forms. */
  validateCreateFromNCBILocation: function(requestData) {
    if (!requestData['refGenomeLabel'].length) {
      alert('Please enter a name for this reference genome.');
      return false;
    }

    if (!requestData['refGenomeAccession'].length) {
      alert('Please specify a genome accession number.');
      return false;
    }

    if (!requestData['importFileFormat'] in ['genbank', 'fasta']) {
      alert('Please select format for file.');
      return false;
    }

    return true;
  },

  /**
   * Creates the S3 uploader if DOM target is present.
   *
   * TODO: Make sure this code still works.
   */
  maybeCreateFineS3Uploader: function() {
    this.uploader = this.$("#uploadDiv").fineUploaderS3({
      debug: true,
      request: {
        endpoint: this.$("#uploadDiv").data("endpoint"),
        accessKey: this.$("#uploadDiv").data("accesskey")
      },
      signature: {
        endpoint: this.$("#uploadDiv").data("signature")
      },
      uploadSuccess: {
        endpoint: this.$("#uploadDiv").data("success")
      },
      objectProperties: {
        key: $.proxy(function(fileId) {
          var filename = this.uploader.fineUploader("getName", fileId);
          return "uploads/" + qq.getUniqueId() + "-" + filename;
        }, this)
      },
      retry: {
        enableAuto: true
      },
      chunking: {
        enabled: true
      },
      deleteFile: {
        endpoint: this.$("#uploadDiv").data("delete"),
        enabled: true,
        forceConfirm: true
      },
      callbacks: {
        onError: function(id, name, reason) {
          alert(reason);
        }
      }
    }).on('complete', $.proxy(function(id, name, response, xhr) {
      var sid = xhr.s3file_id;
      $.post(this.$("#uploadDiv").data("import"), 
            {
              's3file_id': sid,
              'refGenomeLabel': this.$("#uploadDiv").find("#refGenomeLabel").val(),
              'importFileFormat': $("#uploadDiv").find("input[name=importFileFormat]:checked").val(),
            },
            function(data) {
            }
      );}, this));
  },
});
