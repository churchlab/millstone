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

    // For some views, the data is being updated on the server (e.g. alignment
    // is running), so we set this boolean to indicate to clients (i.e.
    // the parent view of this component) that they should consider
    // refreshing with some regular interval while this variable is true.
    // The bit is updated on data fetches.
    this.clientShouldRefresh = false;

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

          if ('clientShouldRefresh' in response) {
            this.clientShouldRefresh = response.clientShouldRefresh;
          }

          this.render();

        }, this));
  },

  getDataDictList: function() {
    // Returns a list in which each element is a dictionary
    // with keys for each of the mData values in fieldConfig
    // and their corresponding value for each object in objList.
    var fields = this.options.fieldConfig;
    var objs = this.options.objList;

    var obj_dict;
    var out_list = [];

    for(var o = 0; o < objs.length; o++) {

      obj_dict = {};

      for(var f = 0; f < fields.length; f++) {

        obj_dict[fields[f].mData] = objs[o][fields[f].mData];
      }

      out_list.push(obj_dict);
    }

    return out_list;
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

    var datatableParams = this.getDataTableParams();

    _.extend(datatableParams, {
      /**********************************************************************
       * Data
       *********************************************************************/
      'aaData': objList,
      'aoColumns': fieldConfig
    });

    if (this.options.extraDatatableParams) {
      _.extend(datatableParams, this.options.extraDatatableParams);
    }

    this.datatable = $('#' + this.datatableId).dataTable(datatableParams);

    // Initialize options for action dropdown menu (next to master checkbox).
    this.createMasterCheckbox();
    this.listenToCheckboxes();
  }
});
