/**
 * @fileoverview List of variants.
 */


gd.AlignmentCreateView = Backbone.View.extend({
  el: '#gd-page-container',


  /** Override. */
  initialize: function() {
    this.render();
  },


  /** Override. */
  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.refGenomeDataTable = new gd.DataTableComponent({
        el: $('#gd-datatable-ref_genome-hook'),
        objList: REF_GENOME_LIST_DATA['obj_list'],
        fieldConfig: REF_GENOME_LIST_DATA['field_config']
    });

    this.samplesDatatable = new gd.DataTableComponent({
        el: $('#gd-datatable-samples-hook'),
        objList: SAMPLES_LIST_DATA['obj_list'],
        fieldConfig: SAMPLES_LIST_DATA['field_config']
    });
  },


  events: {
    'click #gd-align-create-submit-btn': 'handleSubmitClick',
  },


  /**
   * Handles a click on the submit button, aggregating the config
   * options selected and sending a request to the server to start
   * an alignment run.
   */
  handleSubmitClick: function(){
    // Post to this same view for now.
    var postUrl = window.location.pathname;

    // Grab the selected rows.
    var postData = {
        'refGenomes': this.refGenomeDataTable.getCheckedRowUids(),
        'samples': this.samplesDatatable.getCheckedRowUids()
    };

    var onSuccess = function(data) {
      // Redirect according to the server response.
      window.location.href = data.redirect;
    };

    // Execute the post. Should return a redirect response.
    $.post(postUrl, postData, onSuccess, 'json');
  }
});
