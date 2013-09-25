/**
 * @fileoverview Component that takes raw json from the backend intended for
 *     display by jquery.datatable.js and turns the data into a form that
 *     can be rendered as an interactive table.
 */


gd.DataTableComponent = Backbone.View.extend({
  // NOTE: Clients should pass in an element to decorate.

  /** Override. */
  initialize: function() {
    // Handle that will store the reference to the datatable.
    this.datatable = null;

    // Objects converted into displayable form.
    this.displayableObjList = this.makeDisplayableObjectList(
        this.options.objList);

    // Display fields.
    this.displayableFieldConfig = this.makeDisplayableFieldConfig(
        this.options.fieldConfig);

    this.render();
  },


  /** Override. */
  render: function() {
    // Draw the Datatable.
    this.updateDatatable(this.displayableObjList, this.displayableFieldConfig);
  },


  /** Used for updating an already rendered datatable with new data. */
  update: function(newObjList) {
    this.displayableObjList = this.makeDisplayableObjectList(newObjList);
    this.updateDatatable(this.displayableObjList, this.displayableFieldConfig);
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


  /** Make the list of objects into a displayable form. */
  makeDisplayableObjectList: function(objList) {
    var displayableObjList = [];

    _.each(objList, function(obj) {
      var displayableObj = {}

      // Parse through the fields. For any nested objects:
      //     If there is an 'href' and a 'label' key, create an anchor link.
      //     Else if there is only a 'label' key, use only the label.
      //     Else we don't know how to handle this so no change.
      _.each(_.pairs(obj), function(pair) {
        var key = pair[0];
        var value = pair[1];
        var displayValue = value;
        if (_.isArray(value)) {
          displayValue = _.map(value, this.makeDisplayableObject).join(' ');
        }
        else if (typeof(value) == 'object') {
          displayValue = this.makeDisplayableObject(value)
        }
        displayableObj[key] = displayValue;
      }, this);

      // For the object itself, if there is a label and href key,
      // change value keyed by value to be an anchor.
      if ('label' in obj && 'href' in obj) {
        displayableObj['label'] =
            '<a href="' + obj.href + '">' + obj.label + '</>';
      }

      // Add a checkbox.
      var uid = 'undefined';
      if ('uid' in displayableObj) {
        uid = displayableObj['uid'];
      }

      displayableObj['checkbox'] =
          '<input type="checkbox" class="gd-dt-cb" name="gd-row-select" value="' + uid + '">';

      displayableObjList.push(displayableObj);
    }, this);

    return displayableObjList;
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
        'sTitle': 'Select',
        'sClass': 'gd-dt-cb-div',
        'sWidth': '10%',
    });

    return displayableFieldConfig;
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


  /**
   * Finds the gd-dt-master-cb checkbox class and draws a master checkbox
   * and dropdown button that can toggles all checkboxes that are in its table
   */
  createMasterCheckbox: function() {
    $(".gd-dt-cb.master").empty();
    $(".gd-dt-cb.master").append(
      '<div class="gd-dt-cb-div master pull-right btn-group">' +
        '<button class="btn"><input type="checkbox" class="gd-dt-cb master" id="' + this.datatableId + '-master-cb"></button>' +
        '<button class="btn dropdown-toggle" style="min-height: 26px" data-toggle="dropdown">' +
          '<span><i class="icon-chevron-down"></i></span>' +
        '</button>' +
        '<ul class="dropdown-menu" id="' + this.datatableId + '-dropdown">' +
        '</ul>' +
      '</div>');

    this.master_cb = $('#' + this.datatableId + '-master-cb');

    /**
     * If the master checkbox is changed, toggle all checkboxes in the
     * associated table with the following listener.
     */

    this.master_cb.change(_.bind(function(el) {
      // Find all checkboxes in the associated table
      var all_cbs = this.datatable.find('input:checkbox.gd-dt-cb');

      // If none or some of the checkboxes (but not all), then check them all.
      // If all are checked, then uncheck them all.

      var all_checked = _.every(all_cbs, function(cb) {
          return $(cb).is(':checked');})

     _.each(all_cbs, function(cb) {
            $(cb).prop('checked', !all_checked);
            $(cb).triggerHandler('change');
      });

    }, this));
  },


  /**
   * Add a dropdown option to the datatable.
   */
  addDropdownOption: function (html, click_event) {
    var rendered = '<li role="presentation"><a role="menuitem" tabindex="-1" onclick="'+ click_event + '">' + html + '</a></li>';
    $('#' + this.datatableId + '-dropdown').append(rendered);
  },


  /**
   * Updates the datatable view based on the data available.
   *
   * @param {array} objList List of objects to display
   * @param {array} fieldConfig List of column config objects. These must
   *    have the keys (NOTE camel-case):
   *        mData: key corresponding to key in data.
   *        sTitle: title for the column.
   */
  updateDatatable: function(objList, fieldConfig) {
    if (this.datatable == null) {
      // First time we are drawing this table.

      // Create a unique id for the datatable.
      this.datatableId = this.$el.attr('id') + '-datatable';
      this.$el.html(
          '<table cellpadding="0" cellspacing="0" border="0" '+
              'class="table table-striped table-bordered"' +
              'id=' + this.datatableId + '>' +
          '</table>');

      this.datatable = $('#' + this.datatableId).dataTable({
          'aaData': objList,
          'aoColumns': fieldConfig,
          'sDom': "<'row-fluid'<'span3'l><'span5'f><'span3'C><'gd-dt-cb master pull-right span1'>t<'row-fluid'<'span6'i><'span6'p>>",
          'bSortClasses': false,
          'bAutoWidth': false,
          'iDisplayLength': 100,
          'aLengthMenu': [[10, 50, 100, 500, -1], [10, 50, 100, 500, "All"]],
          'sPaginationType': 'bootstrap',
          'fnCreatedRow': this.listenToCheckboxes()
      });
    } else {
      // We are redrawing an existing datatable.
      this.datatable.fnClearTable();
      this.datatable.fnAddData(objList);
    }

    // Initialize options for action dropdown menu (next to master checkbox).
    this.createMasterCheckbox();
    this.listenToCheckboxes();
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
