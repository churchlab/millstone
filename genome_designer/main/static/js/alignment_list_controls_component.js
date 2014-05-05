/**
 * @fileoverview Component that decorates the controls for the list of
 *     ExperimentSamples.
 */


gd.AlignmentListControlsComponent = gd.DataTableControlsComponent.extend({
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
    var deleteOptionHtml = '<a href="#" class="gd-id-samples-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-samples-delete').click(_.bind(this.handleDelete, this));
  },

  /** Show loading feedback while loading. */
  setUIStartLoadingState: function() {
    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin()
  },

  /** Reset UI changes after loading complete.. */
  setUIDoneLoadingState: function() {
    if (this.loadingSpinner) {
      this.loadingSpinner.stop();
    }
  },

  /** Sends request to delete selected samples. */
  handleDelete: function() {
    var uidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!uidList.length) {
      alert("Please select samples to delete.");
      return;
    }

    // Get confirmation from user.
    var agree = confirm("Are you sure you want delete these samples?");
    if (!agree) {
      return;
    }

    this.setUIStartLoadingState();

    var postData = {
        uidList: uidList,
    };

    $.post('/_/alignmentgroups/delete', JSON.stringify(postData),
        _.bind(this.handleDeleteResponse, this));
  },

  handleDeleteResponse: function(response) {
    this.setUIDoneLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
  }
});
