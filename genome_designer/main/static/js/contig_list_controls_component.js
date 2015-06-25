/**
 * @fileoverview Component that decorates the controls for the list of
 *     Contigs.
 */


gd.ContigControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);
    this.decorateControls();
  },

  decorateControls: function() {
    $('#gd-contig-assemble-submit').click(
        _.bind(this.handleAssembleContigs, this));

    this.drawDropdownOptions();
  },

   /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete samples.
    var deleteOptionHtml =
        '<a href="#" class="gd-id-contigs-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-contigs-delete').click(_.bind(this.handleDelete, this));
  },

  /** Send request to generate contigs with default parameters **/
  handleAssembleContigs: function() {
    var requestData = this.prepareRequestData('gd-contig-assemble');
    sample_alignment_uid = requestData['sample_selector']

    var postData = {
        sampleAlignmentUid: sample_alignment_uid
    };

    this.enterLoadingState();

    $.get('/_/alignments/generate_contigs', postData,
        _.bind(this.handleGenerateContigsResponse, this));
  },

  handleGenerateContigsResponse: function(response) {
    this.exitLoadingState();
    $('.modal').modal('hide');

    if (response.is_contig_file_empty == 1) {
      alert('No evidence for structural variants in this alignment');
    } else {
      this.trigger('MODELS_UPDATED');
    };
  },

  /** Parses the form files and prepares the data. */
  prepareRequestData: function(formId) {
    var requestData = {}

    // This works correctly for all inputs except radio buttons
    var formInputs = $('#' + formId + ' :input');
    _.each(formInputs, function(inputObj) {
      requestData[inputObj.name] = inputObj.value;
    });

    return requestData;
  },

  /** Sends request to delete selected samples. */
  handleDelete: function() {
    var contigUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!contigUidList.length) {
      alert('Please select contigs to delete.');
      return;
    }

    // Get confirmation from user.
    var agree = confirm('Are you sure you want delete these contigs?')
    if (!agree) {
      return;
    }

    this.enterLoadingState();

    var postData = {
        contigUidList: contigUidList,
    };

    $.post('/_/contigs/delete', JSON.stringify(postData),
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

  destroy: function() {
    $('#gd-contig-assemble-submit').unbind('click');
  }
});
