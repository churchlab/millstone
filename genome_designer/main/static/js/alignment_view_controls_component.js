/**
 * @fileoverview Component that decorates the controls for the view of a
 *     specific AlignmentGroup and the aligned ExperimentSamples.
 */

gd.AlignmentViewControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);

    this.render();
  },

  render: function() {
    this.decorateControls();
  },

  decorateControls: function() {
    this.drawDropdownOptions();
  },

  /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete samples.
    var downloadBamOptionHTML =
        '<a href="#" class="gd-id-ag-download-bam">Download BAM</a>';
    this.addDropdownOption(downloadBamOptionHTML);
    $('.gd-id-ag-download-bam').click(_.bind(this.handleDownloadBam, this));
  },

  handleDownloadBam: function() {
    var uidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (uidList.length != 1) {
      alert("Please select one sample at a time.");
      return;
    }

    var formJqueryObj = $('#gd-ag-download-bam-form');

    // Reset the form html
    formJqueryObj.empty();

    // Append the form fields.
    this._appendInputFieldToForm(formJqueryObj, 'estaUid', uidList[0]);

    // Submit the form. This cause a download to start.
    formJqueryObj.submit();
  },

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },
});
