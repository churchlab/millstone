/**
 * @fileoverview Component that allows the user to change which fields are
 *     displayed in the Variants DataTable view.
 */


gd.VariantFieldSelectComponent = Backbone.View.extend({
  /**
   * Backbone sugar for specificying element which this component is attached
   * to.
   */
  el: '#gd-filter-key-modal',


  /** Override. */
  initialize: function() {
    this.render();
  },


  /** Override. */
  render: function() {
    $('#gd-filter-key-modal').modal();

    // Render all the key-value pairs.
    this.renderKeyMapSection(this.model.get('variantKeyMap').snp_caller_common_data);
    this.renderKeyMapSection(this.model.get('variantKeyMap').snp_evidence_data);

    $('#gd-filter-key-modal-dismiss-btn').focus();
  },


  renderKeyMapSection: function(key_map) {
    _.each(_.pairs(key_map), function(pair) {
        var maybeChecked = '';
        if ('checked' in pair[1]) {
          maybeChecked = ' checked';
        }

        var rowHtmlString =
            '<tr>' +
              '<td>' + pair[0] + '</td>' +
              '<td>' +
                '<input class="gd-id-filter-key-modal-checkbox" ' +
                    'type="checkbox" value="' + pair[0]+ '"' +
                    maybeChecked + '>' +
              '</td>' +
            '</tr>';
        $('#gd-filter-key-modal-tbody').append(rowHtmlString);
    });
  },


  /** Backbone sugar for registering event listeners. */
  events: {
    'click #gd-filter-key-modal-apply-btn': 'handleApplySelect',
  },


  /** Hides the modal. */
  hide: function() {
    $('#gd-filter-key-modal').modal('hide');
  },


  /**
   * Handles the user indicating that the new field selections should be
   * applied.
   */
  handleApplySelect: function() {
    // Gather the fields that are checked.
    this.trigger('updateSelectedFields', this.getCheckedRowKeyNames());
  },


  /** Returns an array of keys that are selected. */
  getCheckedRowKeyNames: function() {
    var selectedKeyNames = [];
    _.each($('.gd-id-filter-key-modal-checkbox'), function(checkboxEl) {
      if (checkboxEl.checked) {
        selectedKeyNames.push(checkboxEl.value);
      }
    });
    return selectedKeyNames;
  }
});
