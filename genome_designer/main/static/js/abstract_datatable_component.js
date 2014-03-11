/**
 * @fileoverview Abstract base class that wraps the jQuery DataTables plugin.
 *     Clients should use one of the components that inherits from this:
 *         * DataTableComponent: Standard component that has all its data.
 *         * ServerSideDataTableComponent: Pagination handled by server calls.
 */

gd.AbstractDataTableComponent = Backbone.View.extend({

  /** Destroys the datable. */
  destroy: function() {
    if (this.datatable) {
      this.datatable.fnDestroy(true);
    }
    this.datatable = null;
  },

  /** Returns a string for the displayable representation of obj. */
  makeDisplayableObject: function(obj) {
    /* Compute href for object with class information. */
    if ('href' in obj && 'label' in obj) {

      /* If displayable object has css classes, compile them into string */
      if ('classes' in obj) {
        class_str = 'class="' + obj.classes.join(' ') + '" ';
      }
      else {
        class_str = '';
      }

      displayValue = '<a ' + class_str + 'href="' + obj.href + '">' + obj.label + '</>';
    } else if ('label' in obj) {
      displayValue = obj.label;
    } else {
      displayValue = obj;
    }
    return displayValue
  },

  /** Make the field config displayable. */
  makeDisplayableFieldConfig: function(fieldConfig) {
    var displayableFieldConfig = [];

    // Perform a deep copy.
    _.each(fieldConfig, function(col) {
      displayableFieldConfig.push(_.clone(col));
    })

    // Add a column for a select checkbox.
    displayableFieldConfig.push({
        'mData': 'checkbox',
        'sTitle': ' ',
        'sClass': 'gd-dt-cb-div',
        'sWidth': '30px',
    });

    return displayableFieldConfig;
  },

  /**
   * Add a dropdown option to the datatable.
   */
  addDropdownOption: function (html, click_event) {
    var rendered = '<li role="presentation"><a role="menuitem" tabindex="-1" onclick="'+ click_event + '">' + html + '</a></li>';
    $('#' + this.datatableId + '-dropdown').append(rendered);
  },

  /**
   * Listen to newly made datatables checkboxes to update class info and
   * make their parent td clickable.
   */
  listenToCheckboxes: function() {
    $('td.gd-dt-cb-div').click(function(ev){
      if (ev.target.nodeName === "INPUT"){
        return;
      }

      var cb = $(this).find('input:checkbox.gd-dt-cb').get(0);

      if ($(cb).is(':checked')){
        $(cb).prop('checked', false);
      } else {
        $(cb).prop('checked', true);
      }
      $(cb).triggerHandler('change');
    });

    $('input:checkbox.gd-dt-cb').change(function() {
      if ($(this).is(':checked')) {
        $(this).parent('td').addClass('active');
      } else {
        $(this).parent('td').removeClass('active');
      }
    });
  },

  /** Returns an array of uids for the rows that are selected. */
  getCheckedRowUids: function() {
    var selectedUids = [];
    _.each($('input', this.datatable.fnGetNodes()), function(checkboxEl) {
      if (checkboxEl.checked) {
        selectedUids.push(checkboxEl.value);
      }
    });
    return selectedUids;
  }
});
