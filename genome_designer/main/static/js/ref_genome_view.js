/**
 * @fileoverview Single ReferenceGenome view.
 */


gd.RefGenomeView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-refgenomes').addClass('active');
    $('#gd-ref-genome-download-fasta').click(
        _.bind(this.handleDownload, this, 'fasta'));
    $('#gd-ref-genome-download-genbank').click(
        _.bind(this.handleDownload, this, 'genbank'));
    this.redrawDatatable();
  },

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },

  /** Starts a download of the reference genome in the file format passed. */
  handleDownload: function(fileFormat) {
    var formJqueryObj = $('#gd-download-form');

    // Reset the form html
    formJqueryObj.empty();

    // Append the form fields.
    this._appendInputFieldToForm(formJqueryObj,
        'file_format', fileFormat);
    this._appendInputFieldToForm(formJqueryObj,
        'reference_genome_uid', reference_genome_uid);

    // Submit the form. This cause a download to start.
    formJqueryObj.submit();
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    var requestData = {referenceGenomeUid: this.model.get('uid')};

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-ref-genome-view-datatable-hook'),
        serverTarget: '/_/single_ref_genome',
        requestData: requestData,
    });
  }
});
