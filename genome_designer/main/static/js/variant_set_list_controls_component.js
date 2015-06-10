/**
 * @fileoverview Component that decorates the controls for the list of
 *     Variant Sets.
 */


gd.VariantsetsControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);

    if (!this.model) {
      throw "VariantSetControlsComponent requires model.";
    }
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
    var deleteOptionHtml = (
        '<a href="#" class="gd-id-variant-sets-delete">Delete</a>');
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-variant-sets-delete').click(_.bind(this.handleDelete, this));
  },

  /** Sends request to delete selected samples. */
  handleDelete: function() {
    var variantSetUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!variantSetUidList.length) {
      alert("Please select variant sets to delete.");
      return;
    }

    // Get confirmation from user.
    var agree = confirm("Are you sure you want delete these variant sets?");
    if (!agree) {
      return;
    }

    this.setUIStartLoadingState();

    var postData = {
        variantSetUidList: variantSetUidList,
    };

    $.post('/_/variants/delete', JSON.stringify(postData),
        _.bind(this.handleDeleteResponse, this));
  },

  handleDeleteResponse: function(response) {
    this.setUIDoneLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
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

  destroy: function() {
    $('#gd-variant-set-form-from-file-submit').unbind('click')
    $('#gd-variant-set-form-empty-submit').unbind('click')
    this.remove();
    this.unbind();
  }
});
