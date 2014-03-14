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

    if (this.options.hasOwnProperty('serverTarget')) {
      this.initializeFromServerTarget();
    } else {

      this.displayableObjList = this.makeDisplayableObjectList(
          this.options.objList);

      // Display fields.
      this.displayableFieldConfig = this.makeDisplayableFieldConfig(
          this.options.fieldConfig);

      this.render();
    };

  },

  initializeFromServerTarget: function() {
    // If passed a serverTarget, then get the ObjList and FieldConfig
    // from the url via JSON.
    $.get(this.options.serverTarget, this.options.requestData, 
        _.bind(function(response) {
          this.options.objList = response.obj_list;
          this.options.fieldConfig = response.field_config;

          this.displayableObjList = this.makeDisplayableObjectList(
              this.options.objList);

          this.displayableFieldConfig = this.makeDisplayableFieldConfig(
              this.options.fieldConfig);

          this.render();

        }, this));
  },

  /** Override. */
  render: function() {
    // Draw the Datatable.
    this.updateDatatable(this.displayableObjList, this.displayableFieldConfig);

    // Draw controls if the template was supplied.
    if (this.options.hasOwnProperty('controlsTemplate')) {
      this.addControlsFromTemplate(this.options.controlsTemplate, this.options.requestData);
    };

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
    this.$el.append(
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
          "<'gd-datatable-top-container'" +    // first make the top container
          "<'gd-datatable-controlbox'>" +  // model-specific control buttons
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
