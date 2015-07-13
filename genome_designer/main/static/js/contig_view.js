/**
 * @fileoverview Single ReferenceGenome view.
 */


gd.ContigView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-contig-download-fasta').click(
        _.bind(this.handleDownload, this));
  },

  /** Starts a download of the reference genome in the file format passed. */
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

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },

});
