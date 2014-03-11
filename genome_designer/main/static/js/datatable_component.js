/**
 * @fileoverview Component that takes raw json from the backend intended for
 *     display by jquery.datatable.js and turns the data into a form that
 *     can be rendered as an interactive table.
 */


gd.DataTableComponent = gd.AbstractDataTableComponent.extend({
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

    // Activate selectpicker on the dropdowns.
    length_select = $("div.dataTables_length > label > select")
    length_select.addClass('selectpicker')
    length_select.selectpicker({
      width: 'auto'
    });
  },


  /** Used for updating an already rendered datatable with new data. */
  update: function(newObjList, newFieldConfig) {
    this.displayableObjList = this.makeDisplayableObjectList(newObjList);
    this.displayableFieldConfig = this.makeDisplayableFieldConfig(newFieldConfig);
    this.render();
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


  /**
   * Finds the gd-dt-master-cb checkbox class and draws a master checkbox
   * and dropdown button that can toggles all checkboxes that are in its table
   */
  createMasterCheckbox: function() {
    var masterCheckboxElId = this.datatableId + '-master-cb';
    this.$el.find(".gd-dt-cb.master").empty();
    this.$el.find(".gd-dt-cb.master").append(
      '<div class="gd-dt-cb-div master pull-right btn-group">' +
        '<button class="btn btn-default"><input type="checkbox" class="gd-dt-cb master" id="blah"></button>' +
        '<button class="btn btn-default dropdown-toggle" style="min-height: 26px" data-toggle="dropdown">' +
          '<span><i class="caret"></i></span>' +
        '</button>' +
        '<ul class="dropdown-menu" id="' + this.datatableId + '-dropdown">' +
        '</ul>' +
      '</div>');

    /**
     * If the master checkbox is changed, toggle all checkboxes in the
     * associated table with the following listener.
     */

    $('#' + masterCheckboxElId).change(_.bind(function(el) {
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
   * Updates the datatable view based on the data available.
   *
   * @param {array} objList List of objects to display
   * @param {array} fieldConfig List of column config objects. These must
   *    have the keys (NOTE camel-case):
   *        mData: key corresponding to key in data.
   *        sTitle: title for the column.
   */
  updateDatatable: function(objList, fieldConfig) {
    // Clear the existing dattable.
    if (this.datatable != null) {
      this.datatable.fnClearTable();
    }

    // Draw the table.
    // Create a unique id for the datatable.
    this.datatableId = this.$el.attr('id') + '-datatable';
    this.$el.html(
        '<table cellpadding="0" cellspacing="0" border="0" '+
            'class="table table-striped table-bordered"' +
            'id=' + this.datatableId + '>' +
        '</table>');

    var datatableParams = {
        'aaData': objList,
        'aoColumns': fieldConfig,
        'sDom': 
          // This disgusting string tells where to put the table subelements.
          // http://datatables.net/ref#sDom
          // l is the row listing
          // C is the ColVis plugin
          // t is the table
          // i is the row info 
          // p is pagination
          // gd-dt-cb is our master checkbox
          "<'panel panel-default'" +      // start panel containing table
          "<'panel-body'" +               // start panel containing table
          "<'gd-table-controls'" +    // first make the top container
          "<'gd-new-button'>" +  // green button
          "l" +               // records per row
          "<'gd-dt-cb master pull-right'>" + // master cb
          ">>" +  // close panel body, container, row
            "t" + // THE TABLE
          "<'panel-footer'<'container-fluid'<'row'" +    // bottom container/row
            "<'col-sm-4'i><'col-sm-8'p>" + // info, pagination
          ">>>>", // close row, bottom container, panel footer, panel

        'bFilter': false,
        'bSortClasses': false,
        'bAutoWidth': false,
        'iDisplayLength': 100,
        'aLengthMenu': [[10, 50, 100, 500, -1], [10, 50, 100, 500, "All"]],
        'sPaginationType': 'bootstrap',
        'fnCreatedRow': this.listenToCheckboxes()
    };

    if (this.options.extraDatatableParams) {
      _.extend(datatableParams, this.options.extraDatatableParams);
    }

    this.datatable = $('#' + this.datatableId).dataTable(datatableParams);

    // Initialize options for action dropdown menu (next to master checkbox).
    this.createMasterCheckbox();
    this.listenToCheckboxes();
  }
});
