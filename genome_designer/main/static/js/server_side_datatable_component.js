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


gd.ServerSideDataTableComponent = Backbone.View.extend({
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


  /** Destroys the datable. */
  destroy: function() {
    if (this.datatable) {
      this.datatable.fnDestroy(true);
    }
    this.datatable = null;
  },


  /** Used for updating an already rendered datatable with new data. */
  update: function(newObjList, newFieldConfig, numTotalVariants) {
    this.displayableObjList = this.makeDisplayableObjectList(newObjList);
    this.displayableFieldConfig = this.makeDisplayableFieldConfig(
        newFieldConfig);
    this.numTotalVariants = numTotalVariants;
    this.render();
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
    if (!('uid' in displayableObj)) {
      return 'undefined';
    }

    var value = displayableObj.uid;
    if ('sample_uid' in displayableObj) {
      value += ',' + displayableObj.sample_uid;
    }
    return value;
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
    this.$el.find(".gd-dt-cb.master").empty();
    this.$el.find(".gd-dt-cb.master").append(
      '<div class="gd-dt-cb-div master pull-right btn-group">' +
        '<button class="btn"><input type="checkbox" class="gd-dt-cb master" id="' + this.datatableId + '-master-cb"></button>' +
        '<button class="btn dropdown-toggle" style="min-height: 26px" data-toggle="dropdown">' +
          '<span><i class="icon-chevron-down"></i></span>' +
        '</button>' +
        '<ul class="dropdown-menu" id="' + this.datatableId + '-dropdown">' +
        '</ul>' +
      '</div>');

    this.master_cb = $('#' + this.datatableId + '-master-cb');

    // If the master checkbox is changed, toggle all checkboxes in the
    // associated table with the following listener.
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

    var datatableParams = {
        /**********************************************************************
         * Data
         *********************************************************************/
        'aaData': objList,
        'aoColumns': fieldConfig,

        /**********************************************************************
         * Display
         *********************************************************************/
        // Custom positioning for DataTable control elements.
        'sDom': "<'row-fluid'<'span3'l><'span5'f><'span3'C><'gd-dt-cb master pull-right span1'>t<'row-fluid'<'span6'i><'span6'p>>",
        // Don't show the client-side filter box.
        'bFilter': false,
        // Don't add visual highlights to sorted classes.
        'bSortClasses': false,
        // Don't automatically calculate optimal table and col widths.
        'bAutoWidth': false,

        /**********************************************************************
         * Pagination
         *********************************************************************/
        'sPaginationType': 'bootstrap',
        'bServerSide': true,
        'sAjaxSource': this.serverTarget,
        // Inject custom params to filter variants.
        'fnServerParams': this.serverParamsInjector,
        // Override the server function.
        'fnServerData': _.bind(this.customServerDataFn, this),
        // Locked pagination size for now.
        'bLengthChange': false,
        'iDisplayLength': 100,
        // Defer initial load.
        'iDeferLoading': numTotalVariants,

        /**********************************************************************
         * Misc
         *********************************************************************/
        // Called each time a new row is created.
        'fnCreatedRow': this.listenToCheckboxes()
    };

    if (this.options.extraDatatableParams) {
      _.extend(datatableParams, this.options.extraDatatableParams);
    }

    this.datatable = $('#' + this.datatableId).dataTable(datatableParams);

    // // Initialize options for action dropdown menu (next to master checkbox).
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
        fnCallback(this.cleanServerResponse(json, oSettings));

        this.trigger('DONE_LOADING');
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
