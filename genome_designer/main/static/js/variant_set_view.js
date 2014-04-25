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
    $('#gd-variant-set-print-mage-oligos-modal').on('shown.bs.modal',
        _.bind(function (e) {
          if (!this.printMageOligosModal) {
            this.printMageOligosModal=
                new gd.PrintMageOligosModal({model: this.model});
          }
        }, this));
  },

  events: {
    'click #gd-variant-set-view-export-as-csv': 'handleExportCsv',
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

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },
});
