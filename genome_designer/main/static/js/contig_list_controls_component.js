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
