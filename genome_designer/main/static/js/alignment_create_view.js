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
  handleSubmitClick: function() {
    // Post to this same view for now.
    var postUrl = window.location.pathname;

    // Grab the selected rows.
    var postData = {
        refGenomeUidList: this.refGenomeDataTable.getCheckedRowUids(),
        sampleUidList: this.samplesDatatable.getCheckedRowUids()
    };

    var validationResult = this.validatePostData(postData);
    if (!validationResult.is_success) {
      $('#gd-align-create-submit-error-msg').text(validationResult.error_msg);
      $('#gd-align-create-submit-error').show();
      return;
    }

    var onSuccess = function(data) {
      if ('error' in data) {
        $('#gd-align-create-submit-error-msg').text(data.error);
        $('#gd-align-create-submit-error').show();
        return;
      }

      // Redirect according to the server response.
      window.location.href = data.redirect;
    };

    // Execute the post. Handle the response accordingly.
    $.post(postUrl, JSON.stringify(postData), onSuccess, 'json');
  },


  /**
   * Validate the post data before submitting.
   * @param {object} postData The creation data to be posted to the server.
   * @return {object} Response with keys:
   *     * is_success {boolean} Whether validation succeded.
   *     * error_msg {string} Human-readable description of the error.
   */
  validatePostData: function(postData) {
    if (!postData.refGenomeUidList.length) {
      return {
          is_success: false,
          error_msg: 'Please select at least one reference genome.'
      };
    }

    if (!postData.sampleUidList.length) {
      return {
          is_success: false,
          error_msg: 'Please select at least one sample.'
      };
    }

    return {
          is_success: true,
          error_msg: ''
    };
  }
});
