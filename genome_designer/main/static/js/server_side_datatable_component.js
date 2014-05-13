/**
 * @fileoverview Wrapper for the DataTables jQuery plugin whose primary purpose
 *     is delegating pagination-handling to the server side. This is important
 *     for our application because we may be dealing with Variants that number
 *     in the tens and maybe hundreds of thousands so we can't reasonably
 *     send all this data to the frontend.
 *
 * NOTE: This component was created by copying datatable_component.js at commit
 *     7a070f0b7873d96b2cb5e4c54fb34b3c38d45afb and contains duplicated code at
 *     the moment. These components will likely diverge further.
 */


gd.ServerSideDataTableComponent = gd.AbstractDataTableComponent.extend({
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

    // The server target for updates.
    this.serverTarget = this.options.serverTarget;

    // Function that injects additional params for the server-side request.
    this.serverParamsInjector = this.options.serverParamsInjector;

    this.render();
  },


  /** Override. */
  render: function() {
    // Draw the Datatable.
    this.redrawDatable(this.displayableObjList, this.displayableFieldConfig,
        this.numTotalVariants);
  },


  /** Used for updating an already rendered datatable with new data. */
  update: function(newObjList, newFieldConfig, numTotalVariants) {
    this.displayableObjList = this.makeDisplayableObjectList(newObjList);
    this.displayableFieldConfig = this.makeDisplayableFieldConfig(
        newFieldConfig);
    this.numTotalVariants = numTotalVariants;
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
      var value = this.getCheckboxValueFromDisplayableObj(displayableObj);
      displayableObj['checkbox'] =
          '<input type="checkbox" class="gd-dt-cb" name="gd-row-select"' +
              'value="' + value + '">';

      displayableObjList.push(displayableObj);
    }, this);

    return displayableObjList;
  },


  /**
   * Helper method for preparing the checkbox value to be sent to the
   * server.
   *
   * HACK: Hard-coded for the Variant view.
   */
  getCheckboxValueFromDisplayableObj: function(displayableObj) {
    if (!('UID' in displayableObj)) {
      throw "Unexpected format for displayable object.";
    }

    var value = displayableObj.UID;
    if ('EXPERIMENT_SAMPLE_UID' in displayableObj) {
      value += ',' + displayableObj.EXPERIMENT_SAMPLE_UID;
    }
    return value;
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
  redrawDatable: function(objList, fieldConfig, numTotalVariants) {
    // Clear the existing dattable.
    this.destroy();

    // Draw the table.
    // Create a unique id for the datatable.
    this.datatableId = this.$el.attr('id') + '-datatable';
    this.$el.html(
        '<table cellpadding="0" cellspacing="0" border="0" '+
            'class="table table-striped table-bordered"' +
            'id=' + this.datatableId + '>' +
        '</table>');

    var datatableParams =this.getDataTableParams();

    var serverSideDatatableParams = {
      /**********************************************************************
       * Data
       *********************************************************************/
      'aaData': objList,
      'aoColumns': fieldConfig,
      'bFilter': false,
      // Don't add visual highlights to sorted classes.
      /**********************************************************************
       * Server-side options
       *********************************************************************/
      'bServerSide': true,
      'sAjaxSource': this.serverTarget,
      // Inject custom params to filter variants.
      'fnServerParams': this.serverParamsInjector,
      // Override the server function.
      'fnServerData': _.bind(this.customServerDataFn, this),
      // Locked pagination size for now.
      'bLengthChange': false,
      // Defer initial load.
      'iDeferLoading': numTotalVariants
    };

    // Combine common and server side params
    _.extend(datatableParams, serverSideDatatableParams);

    // Add external params that are passed in
    if (this.options.extraDatatableParams) {
      _.extend(datatableParams, this.options.extraDatatableParams);
    }
    
    this.datatable = $('#' + this.datatableId).dataTable(datatableParams);

    // Draw the entity-specific controls that listen for this trigger
    this.trigger('DONE_TABLE_REDRAW');

    // Draw the control buttons at the top-left of the dable component,
    // if passed, the proper template. 
    if (this.options.controlsTemplate) {
      this.addControlsFromTemplate(this.options.controlsTemplate, 
          this.options.requestData);
    }

    // Initialize options for action dropdown menu (next to master checkbox).
    this.createMasterCheckbox();
    this.listenToCheckboxes();
    
  },


  /**
   * Custom function that the DataTable will call to get data from the
   * server. This is mostly copied from the jquery.datatables.js source
   * but modified at the callback point to prepare the data for display.
   */
  customServerDataFn: function(sUrl, aoData, fnCallback, oSettings) {
    this.trigger('START_LOADING');

    oSettings.jqXHR = $.ajax({
      'url':  sUrl,
      'data': aoData,
      'success': _.bind(function (json) {
        if ( json.sError ) {
          oSettings.oApi._fnLog( oSettings, 0, json.sError );
        }
        $(oSettings.oInstance).trigger('xhr', [oSettings, json]);

        // This is the change to the code copied from jquery.datatable.js.
        // We convert the server json response to what is expected by
        // the rest of the DataTables machinery.
        if (json.error) {
          // Send trivial data on to rest of DataTable.
          fnCallback({aaData: []});
          this.trigger('DATA_FETCH_ERROR', json.error);
        } else {
          fnCallback(this.cleanServerResponse(json, oSettings));
          this.trigger('DONE_LOADING');
        }
      }, this),
      'dataType': "json",
      'cache': false,
      'type': oSettings.sServerMethod,
      'error': function(xhr, error, thrown) {
        if (error == "parsererror") {
          oSettings.oApi._fnLog(oSettings, 0,
              "DataTables warning: JSON data from " +
                  "server could not be parsed. This is caused by a " +
                  "JSON formatting error.");
        }
      }
    });
  },


  /**
   * Takes the server response json and creates an object that Datatables
   * can use to redraw itself.
   *
   * For expected format, see: http://datatables.net/usage/server-side
   */
  cleanServerResponse: function(json, oSettings) {
    var aaData = [];

    // Parse the variant data.
    var variantListData = JSON.parse(json.variant_list_json);
    if (variantListData.obj_list.length) {
      aaData = this.makeDisplayableObjectList(variantListData.obj_list);
    }

    datatablesJson = {
        iTotalRecords: Number(json.num_total_variants),
        iTotalDisplayRecords: Number(json.num_total_variants),
        aaData: aaData,
    };
    return datatablesJson;
  }
});
