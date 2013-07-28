/**
 * @fileoverview Component that takes raw json from the backend intended for
 *     display by jquery.datatable.js and turns the data into a form that
 *     can be rendered.
 */


gd.DataTableComponent = Backbone.View.extend({
  // NOTE: Clients should pass in an element to decorate.

  /** Override. */
  initialize: function() {
    // Raw list of objects.
    this.objList = this.options.objList;

    // Objects converted into displayable form.
    this.displayableObjList = this.makeDisplayableObjectList(this.objList);

    // List of fields.
    this.fieldConfig = this.options.fieldConfig;

    this.render();
  },


  /** Override. */
  render: function() {
    // Draw the Datatable.
    this.updateDatatable(this.displayableObjList, this.fieldConfig);
  },


  /** Make the list of objects into a displayable form. */
  makeDisplayableObjectList: function(objList) {
    var displayableObjList = [];

    _.each(objList, function(obj) {
      displayableObj = {}
      _.each(_.pairs(obj), function(pair) {
        var key = pair[0];
        var value = pair[1];
        var displayValue = value;
        if (typeof(value) == 'object') {
          displayValue = '<a href="' + value.href + '">' + value.label + '</>';
        }
        displayableObj[key] = displayValue;
      })

      displayableObjList.push(displayableObj);
    });

    return displayableObjList;
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
    // Create a unique id for the datatable.
    var datatableId = this.$el.attr('id') + '-datatable';
    this.$el.html(
        '<table cellpadding="0" cellspacing="0" border="0" '+
            'class="table table-striped table-bordered"' +
            'id=' + datatableId + '>' +
        '</table>');
    $('#' + datatableId).dataTable({
        'aaData': objList,
        'aoColumns': fieldConfig,
        "bSortClasses": false,
        'sPaginationType': 'bootstrap'
    });
  }
});
