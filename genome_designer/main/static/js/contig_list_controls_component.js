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
    // Option to delete contigs.
    var deleteOptionHtml =
        '<a href="#" class="gd-id-contigs-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-contigs-delete').click(_.bind(this.handleDelete, this));

    // Option to find contig insertion location.
    var findInsertionLocationOptionHtml =('<a href="#" class="gd-id-' +
        'contigs-find-insertion-location">Find Insertion Location</a>');
    this.addDropdownOption(findInsertionLocationOptionHtml);
    $('.gd-id-contigs-find-insertion-location').click(_.bind(
        this.handleFindInsertionLocation, this));

    // Option to place contig.
    var placeInRefOptionHtml = ('<a href="#" class="gd-id-contigs-place' +
        '-in-ref">Make New Reference with Contig Incorported</a>');
    this.addDropdownOption(placeInRefOptionHtml);
    $('.gd-id-contigs-place-in-ref').click(_.bind(
        this.handlePlaceInRef, this));
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

  /** Sends request to find insertion location. */
  handleFindInsertionLocation: function() {
    var contigUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!contigUidList.length) {
      alert('Please select contigs for which to find insertion locations.');
      return;
    }

    // Get confirmation from user.
    var getData = {
        contigUidList: contigUidList
    };

    $.get('/_/contigs/has_insertion_location', {data:JSON.stringify(getData)},
        _.bind(function(response){
          var hasInsertionEndpoints = response.has_insertion_location;
          if (hasInsertionEndpoints) {
            var agree = confirm('Some of the contigs selected already have' +
                'insertion locations; are you sure you want to overwrite ' +
                'them?')
            if (!agree) {
              return;
            }
          }
          this.enterLoadingState();
          var postData = {
              contigUidList: contigUidList,
          }
          $.post('/_/contigs/find_insertion_location',
              JSON.stringify(postData),
              _.bind(this.handleFindInsertionLocationResponse, this));
        }, this))
  },

  handleFindInsertionLocationResponse: function(response) {
    this.exitLoadingState();
    if ('error' in response) {
      var errorMessage = ('Unable to find insertion locations for the ' +
          'following contigs:\n')
      for (i in response.error) {
        errorMessage += (response.error[i][0] + ': ' + response.error[i][1] +
            '\n')
      }
      var agree = confirm(errorMessage)
    }
    
    this.trigger('MODELS_UPDATED');
  },

  handlePlaceInRef: function() {

    var contigUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!contigUidList.length) {
      alert('Please select a contig to incorporate');
      return;
    }

    // If multiple selected, show message
    if (contigUidList.length > 1) {
      alert('Please select only one contig to incorporate');
      return;
    }

    var getData = {
        contigUidList: contigUidList
    };

    $.get('/_/contigs/has_insertion_location', {data:JSON.stringify(getData)},
        _.bind(function(response){
          var hasInsertionEndpoints = response.has_insertion_location;
          if (!hasInsertionEndpoints) {
            alert('The contig you selected does not yet have ' +
                'insertion location data associated.  Before incorporating ' +
                'it in a new genome, please find insertion location data ' +
                'for it by selecting "Find Insertion Location" from the ' +
                'dropdown')
            return;
          }
          this.handlePlaceInRefLabelPrompt();
        }, this))
  },

  handlePlaceInRefLabelPrompt: function() {
    // Get new genome name
    var newGenomeLabel = prompt(
        'Enter a name for the new genome:', 'new_genome_name');
    while (newGenomeLabel == '') {
      var newGenomeLabel = prompt(
          'Please enter a non-zero length name for the new genome',
          'new_genome_name');
    }
    if (newGenomeLabel == null) {
      return;
    }

    this.enterLoadingState();

    var contigUidList = this.datatableComponent.getCheckedRowUids();
    var postData = {
        'newGenomeLabel': newGenomeLabel,
        'contigUidList': contigUidList,
    };

    $.post('/_/contigs/place_in_ref', JSON.stringify(postData),
        _.bind(this.handlePlaceInRefResponse, this));
  },

  handlePlaceInRefResponse: function(response) {
    this.exitLoadingState();
  },

  /** Sends request to delete selected contigs. */
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
