/**
 * @fileoverview Single alignment and samples within.
 */


gd.AlignmentView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    this.decorateSidebar();

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-alignment_view-datatable-hook'),
        objList: EXPERIMENT_TO_SAMPLE_DATA['obj_list'],
        fieldConfig: EXPERIMENT_TO_SAMPLE_DATA['field_config'],
        controlsTemplate: '/_/templates/alignment_controls',
        requestData: {
            projectUid: this.model.get('project').uid,
            alignmentGroupUid: this.model.get('alignment_group').uid
        },
    });

    // hook up the button to re-running the alignment as soon
    // as they are drawn.
    this.listenTo(this.datatable, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateAlignmentControls, this));
  },

  handleGenerateContigs: function() {
    var sampleAlignmentUidList = this.datatable.getCheckedRowUids();

    // If nothing to do, show message.
    if (!sampleAlignmentUidList.length) {
      alert('Please select Experiment Samples to use for contig generation.');
      return;
    }

    if(sampleAlignmentUidList.length > 1) {
      alert('Please select only one Experiment Sample to use for contig ' +
          'generation.');
      return;
    }

    var postData = {
        sampleAlignmentUid: sampleAlignmentUidList[0]
    };

    this.enterLoadingState();

    $.get('/_/alignments/generate_contigs', postData,
        _.bind(this.handleGenerateContigsResponse, this));
  },

  handleGenerateContigsResponse: function(response) {
    this.exitLoadingState();
    if (response.is_contig_file_empty == 1) {
      alert('No evidence for structural variants in this alignment');
    };
  },

 /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },

  /** Decorate the side nav bar. */
  decorateSidebar: function() {
    $('#gd-sidenav-link-alignments').addClass('active');
  },

  /** Create component to wrap samples controls, and listen for events. */
  decorateAlignmentControls: function() {
    // Draw 'Generate contigs' button only if the feature is enabled
    if($('#gd-generate-contigs-btn').length) {
      $('#gd-generate-contigs-btn').click(_.bind(
          this.handleGenerateContigs, this));
    }

    // Draw 'Go to variants' button only if AlignmentGroup is COMPLETED.
    if (this.model.get('alignment_group').status == 'COMPLETED') {
      var templateData = {
          projectUid: this.model.get('project').uid,
          alignmentGroupUid: this.model.get('alignment_group').uid,
      };
      $('#gd-ag-controls-toolbar').append(_.template(
        '<a id="gd-ag-go-to-variants-btn" class="btn btn-primary" ' +
            'href="/projects/<%= projectUid %>/analyze/' +
                '<%= alignmentGroupUid %>/variants">' +
          'Go to Variants for this alignment&hellip;' +
        '</a>'
      , templateData));
    }

    // Draw 'Re-run' only if AlignmentGroup is COMPLETED or FAILED.
    if (_.contains(['COMPLETED', 'FAILED'],
        this.model.get('alignment_group').status)) {
      $('#gd-ag-controls-toolbar').append(
        '<button id="gd-alignments-rerun-variants-btn" ' +
            'class="btn btn-primary">' +
          'Re-run&hellip;' +
        '</button>'
      );
      $('#gd-alignments-rerun-variants-btn').click(
          _.bind(this.handleRerunVariantsClick, this));
    }

    $('#gd-ag-controls-toolbar').append(
      '<span id="gd-ag-status-text" class="gd-ag-status-text">' +
        this.model.get('alignment_group').status +
      '</span>');
  },

  /**
   * Tells the server to re-run variant calling for this alignment group.
   * This will also re-run any failed alignments.
   */
  handleRerunVariantsClick: function() {
    // Make sure the user confirms they want to do this.
    var agree = confirm("Are you sure you want to re-run alignment? " +
        "This will delete previously-called Variants.");
    if (!agree) {
      return;
    }

    // Disable re-run button and hide status text. Refresh below.
    $('#gd-alignments-rerun-variants-btn').prop('disabled', true);
    $('#gd-alignments-rerun-variants-btn').text('Re-running...');
    $('#gd-ag-status-text').hide();

    // Disable 'go to variants'. Note that it's an anchor so we disable in a
    // slightly different manner.
    $('#gd-ag-go-to-variants-btn').addClass('disabled');

    // Post to this same view for now.
    var postUrl = window.location.pathname;

    // Nothing to send as this is the only post request this page expects.
    var postData = {};

    var onSuccess = function(data) {
      window.location.reload();
    };

    // // Execute the post. Should return a redirect response.
    $.post(postUrl, JSON.stringify(postData), onSuccess, 'json');
  },

    /** Puts UI in the loading state. */
  enterLoadingState: function() {
    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin();
  },

  /** Puts UI in the loading state. */
  exitLoadingState: function() {
    this.loadingSpinner.stop();
  }
});
