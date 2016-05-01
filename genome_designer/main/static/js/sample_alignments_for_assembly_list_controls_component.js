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
    var callFromContigsOptionHtml =
        '<a href="#" class="gd-id-call-sv-from-contigs">Call SVs using assembled contigs</a>';
    this.addDropdownOption(assembleOptionHtml);
    this.addDropdownOption(callFromContigsOptionHtml);
    $('.gd-id-assemble').click(_.bind(this.handleAssembleContigs, this));
    $('.gd-id-call-sv-from-contigs').click(_.bind(this.handleCallSVFromContigs, this));
  },

  /** Send request to generate contigs with default parameters **/
  handleAssembleContigs: function() {

    sample_alignment_uid_list = this.datatableComponent.getCheckedRowUids()

    var postData = {
        sampleAlignmentUidList: sample_alignment_uid_list
    };

    $.post('/_/alignments/generate_contigs', JSON.stringify(postData),
        _.bind(this.handleAssembleContigsResponse, this));
  },

  /** Send request to call svs from contigs with default parameters **/
  handleCallSVFromContigs: function() {

    sample_alignment_uid_list = this.datatableComponent.getCheckedRowUids()

    var postData = {
        sampleAlignmentUidList: sample_alignment_uid_list
    };

    $.post('/_/alignments/call_svs_from_contigs', JSON.stringify(postData),
        _.bind(this.handleAssembleContigsResponse, this));
  },

  handleAssembleContigsResponse: function(response) {

    if (response.is_contig_file_empty == 1) {
      alert('No evidence for structural variants in this alignment');
    } else {
      this.trigger('MODELS_UPDATED');
    };
  },
});
