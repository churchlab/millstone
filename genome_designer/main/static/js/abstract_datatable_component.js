/**
 * @fileoverview Abstract base class that wraps the jQuery DataTables plugin.
 *     Clients should use one of the components that inherits from this:
 *         * DataTableComponent: Standard component that has all its data.
 *         * ServerSideDataTableComponent: Pagination handled by server calls.
 */

gd.AbstractDataTableComponent = Backbone.View.extend(
{
  /** Prevent initialization. */
  initialize: function() {
    throw "Don't initialize abstract class!";
  },

  /** Destroys the datable. */
  destroy: function() {
    if (this.datatable) {
      this.datatable.fnDestroy(true);
    }
    this.datatable = null;
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
    displayableFieldConfig.unshift({
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
  addDropdownOption: function (html, clickEvent) {
    var rendered = '<li role="presentation"><a role="menuitem" tabindex="-1" onclick="' + 
        clickEvent + '">' + html + '</a></li>';
    $('#' + this.datatableId + '-dropdown').append(rendered);
  },

  /**
   * Finds the gd-dt-master-cb checkbox class and draws a master checkbox
   * and dropdown button that can toggles all checkboxes that are in its table
   */
  createMasterCheckbox: function() {
    var masterCheckboxElId = this.datatableId + '-master-cb';
    this.$el.find(".gd-dt-cb.master").empty();
    this.$el.find(".gd-dt-cb.master").append(
      '<div class="gd-dt-cb-div master pull-left btn-group">' +
        '<button class="btn btn-default">' +
        '<input type="checkbox" class="gd-dt-cb master" id="blah"></button>' +
        '<button class="btn btn-default dropdown-toggle"' +
        ' style="min-height: 26px" data-toggle="dropdown">' +
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
  },

  /**
    * Adds html controls (a dropdown button above the table and accompanying
    * modals. Takes a templateURL that is handled by template_xhrs.py and a
    * requestData object with the required data. Both of these should be
    * passed by whatever asks for the table to be drawn, via this.options.
    */
  addControlsFromTemplate: function(templateURL, requestData) {

    // Add the id of this datatable to the request. That way, we'll
    // know where to add the control box in case there are multiple
    // tables.
    requestData['tableId'] = this.$el.attr('id');

    $.get(templateURL, requestData, _.bind(function(response) {

      // get the location for the buttons
      var controls_div = $("#"+this.$el.attr('id')).find('.gd-datatable-controlbox');

      // remove any old controls
      controls_div.find('.gd-datatable-control').remove();

      // First, put everything at the end of the DOM.
      $('body').append(response);

      // Then, move the .gd-datatable-control stuff into the controlbox.
      $('#' + this.$el.attr('id') +'-control').appendTo(controls_div);

      // Some views need to do stuff to these controls afterwards, so
      // let them know we're done.
      this.trigger('DONE_CONTROLS_REDRAW');
    }, this));
  },

  getDataTableParams: function(dataTableObj) {
    var dataTableParams = {
      /**********************************************************************
       * Display
       *********************************************************************/
      'sDom':
        // Custom positioning for DataTable control elements.
        // This disgusting string tells where to put the table subelements.
        // http://datatables.net/ref#sDom
        // l is the row listing
        // C is the ColVis plugin
        // t is the table
        // i is the row info 
        // p is pagination
        // gd-dt-cb is our master checkbox
        "<'panel panel-default'" +         // start panel containing table
        "<'panel-body'" +                  // start panel containing table
        "<'gd-datatable-controlbox'" +     // first make the top container
        //"l" +                              // records per row (skipping)
        "<'gd-dt-cb master pull-left'>" +  // master cb
        "<'pull-right'ip>" +               // info, pagination
        ">>" +                             // close panel body, container, row
          "t" +                            // THE TABLE
        ">>", // close panel footer, panel
      /**********************************************************************
       * Hide verbose pagination labels
       *********************************************************************/
      "oLanguage": {
        "sInfo": "_START_-_END_ of _TOTAL_",
        "oPaginate": {
          "sPrevious": "",
          "sNext": ""
        }
      },
      /**********************************************************************
       * Data
       *********************************************************************/
      'bFilter': false,
      // Don't add visual highlights to sorted classes.
      'bSortClasses': false,
      // Don't automatically calculate optimal table and col widths.
      'bAutoWidth': false,
      /**********************************************************************
       * Pagination
       *********************************************************************/
      'sPaginationType': 'bootstrap',
      'iDisplayLength': 100,
      /**********************************************************************
       * Misc
       *********************************************************************/
      // Called each time a new row is created.
      'fnCreatedRow': this.listenToCheckboxes()
    };
    return dataTableParams;
  }
});
