/**
 * @fileoverview Component that decorates the controls for the list of
 *     ReferenceGenomes.
 */


gd.RefGenomeControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);

    this.decorateControls();

    this.maybeCreateFineS3Uploader();
  },

  decorateControls: function() {
    $('#gd-ref-genome-upload-through-browser-submit').click(
        _.bind(this.handleUploadThroughBrowser, this));
    $('#gd-ref-genome-upload-from-server-location-submit').click(
        _.bind(this.handleCreateFromServerLocation, this));
    $('#gd-ref-genome-create-from-ncbi-submit').click(
        _.bind(this.handleCreateFromNCBI, this));

    this.drawDropdownOptions();
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

  /** Uploads genome through browser. */
  handleUploadThroughBrowser: function() {
    // Parse the inputs.
    var requestData = this.prepareRequestData(
        'gd-ref-genome-upload-through-browser');

    // Validate the request client-side, even though we just grab
    // the data again into FormData form below.
    if (!this.validateUploadThroughBrowserForm(requestData)) {
      return;
    }

    // Put UI in loading state.
    this.enterLoadingState();

    var onSuccess = _.bind(function(responseData) {
      this.exitLoadingState();

      // Check for error and show in ui. Don't reload the page.
      if (responseData.error.length) {
        alert('Error creating reference genome: ' + responseData.error);
        return;
      }

      this.trigger('MODELS_UPDATED');
      $('.modal').modal('hide');
    }, this);

    var formData = new FormData(
        $('#gd-ref-genome-upload-through-browser')[0]);
    $.ajax({
      url: '/_/ref_genomes/upload_through_browser',
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
    var jqxhr = $.post('/_/ref_genomes/create_from_server_location', requestData,
        _.bind(function(responseData) {
          this.exitLoadingState();

          // Check for error and show in ui. Don't reload the page.
          if (responseData.error.length) {
            alert('Error creating reference genome: ' + responseData.error);
            return;
          }

          // Success, emit an event to be handled by parent.
          this.trigger('MODELS_UPDATED');
          $('.modal').modal('hide');
        }, this));

    jqxhr.fail(_.bind(function() {
      this.exitLoadingState()
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

  validateUploadThroughBrowserForm: function(requestData) {
    if (!requestData['refGenomeLabel'].length) {
      alert('Please enter a name for this reference genome.');
      return false;
    }

    if (!requestData['refGenomeFile'].length) {
      alert('Please specify a file location.');
      return false;
    }

  if (!requestData['importFileFormat'] in ['genbank', 'fasta']) {
      alert('Please select format for file.');
      return false;
    }

    return true;
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
    var jqxhr = $.post('/_/ref_genomes/create_from_ncbi', requestData,
        _.bind(function(responseData) {
          this.exitLoadingState();

          // Check for error and show in ui. Don't reload the page.
          if (responseData.error.length) {
            alert('Error creating reference genome: ' + responseData.error);
            return;
          }

          // Success, emit an event to be handled by parent.
          this.trigger('MODELS_UPDATED');
          $('.modal').modal('hide');
        }, this));

    jqxhr.fail(_.bind(function() {
      this.exitLoadingState()
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

   /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete samples.
    var deleteOptionHtml =
        '<a href="#" class="gd-id-refgenomes-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-refgenomes-delete').click(_.bind(this.handleDelete, this));

    var concatenateOptionHtml = 
        '<a href="#" class="gd-id-refgenomes-concatenate">Concatenate</a>';
    this.addDropdownOption(concatenateOptionHtml);
    $('.gd-id-refgenomes-concatenate').click(
        _.bind(this.handleConcatenate, this));
  },

  /** Sends request to delete selected samples. */
  handleDelete: function() {
    var refGenomeUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!refGenomeUidList.length) {
      alert('Please select reference genomes to delete.');
      return;
    }

    // Get confirmation from user.
    var agree = confirm(
        'Are you sure you want delete these reference genomes? ' +
        'This will also delete any alignments and variants associated ' +
        'with this reference genome.');
    if (!agree) {
      return;
    }

    this.enterLoadingState();

    var postData = {
        refGenomeUidList: refGenomeUidList,
    };

    $.post('/_/ref_genomes/delete', JSON.stringify(postData),
        _.bind(this.handleDeleteResponse, this));
  },

  handleDeleteResponse: function(response) {
    this.exitLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
  },

  handleConcatenate: function() {
    var refGenomeUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!refGenomeUidList.length) {
      alert('Please select reference genomes to concatenate.');
      return;
    }
    // If only one selected, show message.
    if (refGenomeUidList.length == 1) {
      alert('Please select more than one reference genome to concatenate.');
      return;
    }

    // Get new genome name
    var newGenomeLabel = prompt(
        'Enter a name for the concatenated genome:', 'new_genome_name');
    while (newGenomeLabel == '') {
      var newGenomeLabel = prompt(
          'Please enter a non-zero length name for the concatenated genome',
          'new_genome_name');
    }
    if (newGenomeLabel == null) {
      return;
    }

    this.enterLoadingState();

    var postData = {
        'newGenomeLabel': newGenomeLabel,
        'refGenomeUidList': refGenomeUidList,
    };

    $.post('/_/ref_genomes/concatenate', {data:JSON.stringify(postData)},
        _.bind(this.handleConcatenateResponse, this));

},

  handleConcatenateResponse: function(response) {
    this.exitLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
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
              'refGenomeLabel': this.$("#uploadDiv").find(
                  "#refGenomeLabel").val(),
              'importFileFormat': $("#uploadDiv").find(
                  "input[name=importFileFormat]:checked").val(),
            },
            function(data) {
            }
      );}, this));
  },
});
