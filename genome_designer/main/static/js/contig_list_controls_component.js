/**
 * @fileoverview Component that decorates the controls for the list of
 *     Contigs.
 */


gd.ContigControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);
    this.render();
  },

  render: function() {
    this.decorateControls();
  },

  decorateControls: function() {
    this.drawDropdownOptions();
    this.listenToControls();
  },

   /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete contigs.
    var deleteOptionHtml =
        '<a href="#" class="gd-id-contigs-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-contigs-delete').click(_.bind(this.handleDelete, this));
  },

  listenToControls: function() {
    $('#gd-contigs-view-download-btn').click(
        _.bind(this.handleDownloadContigsButtonClick, this));
  },

  handleDownloadContigsButtonClick: function() {
    var formJqueryObj = $('#gd-contigs-export-form');

    // Reset the form html
    formJqueryObj.empty();

    // Append the form fields.
    this._appendInputFieldToForm(formJqueryObj, 'alignment_group_uid',
        this.model.get('alignmentGroupUid'));

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
    gd.Spinner.globalClear();
  },

  destroy: function() {
    $('#gd-contig-assemble-submit').unbind('click');
  }
});
