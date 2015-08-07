/**
 * @fileoverview Single Contig view.
 */


gd.ContigView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-contig-download-fasta').click(
        _.bind(this.handleDownload, this));
    $('#gd-contig-jbrowse').click(
        _.bind(this.handleJbrowse, this));
  },

  /** Starts a download of the contig. */
  handleDownload: function() {
    var formJqueryObj = $('#gd-download-form');

    // Reset the form html
    formJqueryObj.empty();

    // Append the form fields.
    this._appendInputFieldToForm(formJqueryObj,
        'contig_uid', contig_uid);

    // Submit the form. This cause a download to start.
    formJqueryObj.submit();
  },

   /** Makes jbrowse tracks for contig and redirects to jbrowse. */
  handleJbrowse: function() {

    var getData = {
        contigUid: INIT_JS_DATA.entity.uid
    };

    this.enterLoadingState();

    $.get('/_/contigs/make_contig_jbrowse_tracks', getData,
          _.bind(function(response) {
              var jbrowseLink = response.jbrowse_link;
              this.exitLoadingState();
              window.location.href = jbrowseLink;
          }, this));
  },

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
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
