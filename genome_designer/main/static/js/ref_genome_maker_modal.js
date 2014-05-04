/**
 * Modal that allows choosing options for printing MAGE oligos.
 */


gd.RefGenomeMakerModal = Backbone.View.extend({

  initialize: function() {
    this.render();
  },

  render: function() {
    this.listenToControls();
  },

  listenToControls: function() {
    $('#gd-variant-sets-ref-genome-maker-submit').click(
        _.bind(this.handleCreateRefGenome, this));
  },

  handleCreateRefGenome: function() {
    this.enterLoadingState();
    var formData = new FormData($('#gd-variant-set-ref-genome-maker-form')[0]);
    $.ajax({
      url: '/_/sets/generate_new_ref_genome',
      type: 'POST',
      data: formData,
      success: _.bind(this.handleCreateRefGenomeSuccess, this),

      // The following 3 param settings are necessary for properly passing
      // formData. See: http://stackoverflow.com/questions/166221/how-can-i-upload-files-asynchronously-with-jquery
      cache: false,
      contentType: false,
      processData: false
    });
  },

  handleCreateRefGenomeSuccess: function(response) {
    if ('error' in response && response.error.length) {
      this.exitLoadingState();
      alert('ERROR: ' + response.error);
    } else {
      window.location.href = response.redirect;
    }
  },

  /** Puts UI in the loading state. */
  enterLoadingState: function() {
    $("#gd-variant-sets-ref-genome-maker-submit").prop('disabled', true);

    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin();
  },

  /** Exits loading state. */
  exitLoadingState: function() {
    $("#gd-variant-sets-ref-genome-maker-submit").prop('disabled', false);
    this.loadingSpinner.stop();
  },
});
