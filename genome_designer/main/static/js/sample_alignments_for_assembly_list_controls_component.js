/**
 * @fileoverview Component that decorates the controls for the list of
 *     ExperimentSampleToAlignments on the Contig Analyze subview
 */


gd.SampleAssemblyControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);
    this.decorateControls();
  },

  decorateControls: function() {
    this.drawDropdownOptions();
  },

   /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete contigs.
    var assembleOptionHtml =
        '<a href="#" class="gd-id-assemble">Assemble Contigs</a>';
    this.addDropdownOption(assembleOptionHtml);
    $('.gd-id-assemble').click(_.bind(this.handleAssembleContigs, this));
  },

  /** Send request to generate contigs with default parameters **/
  handleAssembleContigs: function() {

    sample_alignment_uid_list = this.datatableComponent.getCheckedRowUids()

    var postData = {
        sampleAlignmentUidList: sample_alignment_uid_list
    };

    this.enterLoadingState();

    $.post('/_/alignments/generate_contigs', JSON.stringify(postData),
        _.bind(this.handleAssembleContigsResponse, this));
  },

  handleAssembleContigsResponse: function(response) {
    this.exitLoadingState();
    $('.modal').modal('hide');

    if (response.is_contig_file_empty == 1) {
      alert('No evidence for structural variants in this alignment');
    } else {
      this.trigger('MODELS_UPDATED');
    };
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
});
