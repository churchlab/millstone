/**
 * @fileoverview Single alignment and samples within.
 */


gd.VariantSetView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.filterString = 'variant_set_uid = ' + this.model.get('uid');
    this.render();
  },

  render: function() {
    this.decorateSidebar();

    this.variantsTableComponent = new gd.VariantsTableComponent({
      model: this.model,
      filterString: this.filterString,
      filterDisabled: true,
    });

    this.listenToControls();
  },

  /** Decorate the side nav bar. */
  decorateSidebar: function() {
    $('#gd-sidenav-link-variant-sets').addClass('active');
  },

  listenToControls: function() {
    // Modal for printing MAGE oligos.
    $('#gd-variant-set-print-mage-oligos-modal').on('shown.bs.modal',
        _.bind(function (e) {
          if (!this.printMageOligosModal) {
            this.printMageOligosModal=
                new gd.PrintMageOligosModal({model: this.model});
          }
        }, this));

    // Modal for generating reference genome.
    $('#gd-variant-set-ref-genome-maker-modal').on('shown.bs.modal',
        _.bind(function (e) {
          if (!this.refGenomeMakerModal) {
            this.refGenomeMakerModal =
                new gd.RefGenomeMakerModal({model: this.model});
          }
        }, this));
  },

  events: {
    'click #gd-variant-set-view-export-as-csv': 'handleExportCsv',
    'click #gd-variant-set-view-delete' : 'handleDeleteSet',
  },

  handleDeleteSet: function() {
    var agree = confirm(
      "Are you sure you want to delete this variant set?\n" +
      "Individual variants will be preserved."
      );
    if (!agree) {
      return;
    }

    this.enterLoadingState()

    var variantSetUidList = [this.model.get('uid')];

    var postData = {
      variantSetUidList: variantSetUidList,
    }
    $.post('/_/variants/delete',JSON.stringify(postData),
      _.bind(this.handleDeleteResponse,this));
  },

  handleDeleteResponse: function(response) {
    this.exitLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      var projectUid = this.model.get('project').uid;
      window.location.href = "/projects/" + projectUid + "/sets";
    }
  },

  /** Dynamically populate form fields and submit. */
  handleExportCsv: function() {
    var formJqueryObj = $('#gd-variant-sets-export-csv-form');

    // Append the form fields.
    this._appendInputFieldToForm(formJqueryObj, 'ref_genome_uid',
        this.model.get('refGenomeUid'));
    this._appendInputFieldToForm(formJqueryObj, 'filter_string',
        this.filterString);

    // Submit the form. This cause a download to start.
    formJqueryObj.submit();
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

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },
});
