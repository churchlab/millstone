/**
 * @fileoverview Component that takes raw json from the backend intended for
 *     display by jquery.datatable.js and turns the data into a form that
 *     can be rendered.
 */


gd.DataTableComponent = Backbone.View.extend({
  // NOTE: Clients should pass in an element to decorate.

  /** Override. */
  initialize: function() {
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


  /** Make the list of objects into a displayable form. */
  makeDisplayableObjectList: function(objList) {
    var displayableObjList = [];

    _.each(objList, function(obj) {
      var displayableObj = {}

      // For any nested objects, try to render them as links.
      _.each(_.pairs(obj), function(pair) {
        var key = pair[0];
        var value = pair[1];
        var displayValue = value;
        if (typeof(value) == 'object') {
          displayValue = '<a href="' + value.href + '">' + value.label + '</>';
        }
        displayableObj[key] = displayValue;
      })

      // For the object itself, if there is a label and href key,
      // change value keyed by value to be an anchor.
      if ('label' in obj && 'href' in obj) {
        displayableObj['label'] =
            '<a href="' + obj.href + '">' + obj.label + '</>';
      }

      // Add a checkbox.
      displayableObj['checkbox'] = '<input type="checkbox">';

      displayableObjList.push(displayableObj);
    });

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
        'sTitle': 'Select'
    });

    return displayableFieldConfig;
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
        'sDom': "<'row'<'span5'l><'span5'f><'align-right span2'C>r>t<'row'<'span6'i><'span6'p>>",
        "bSortClasses": false,
        'sPaginationType': 'bootstrap'
    });
  }
});
