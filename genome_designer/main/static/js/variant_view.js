/**
 * @fileoverview List of variants.
 */


gd.VariantView = Backbone.View.extend({
  el: '#gd-page-container',


  initialize: function() {
    this.render();

    this.updateDatatable();
  },


  render: function() {
    $('#gd-sidenav-link-variants').addClass('active');
  },


  /** Update the datatable. */
  updateDatatable: function() {
    $('#gd-datatable-hook').html(
        '<table cellpadding="0" cellspacing="0" border="0" class="display"' +
            'id="gd-datatable">' +
        '</table>');
    $('#gd-datatable').dataTable( {
        "aaData": [
            {
              "engine": "Trident",
              "browser": "Internet Explorer 4.0",
              "platform": "Win 95+",
              "version": "4",
              "grade": "X"
            },
            {
              "engine": "Trident",
              "browser": "Internet Explorer 5.0",
              "platform": "Win 95+",
              "version": "5",
              "grade": "C"
            },
            {
              "engine": "Trident",
              "browser": "Internet Explorer 5.5",
              "platform": "Win 95+",
              "version": "5.5",
              "grade": "A"
            }
        ],
        "aoColumns": [
            { "mData": "engine", "sTitle": "Engine Blah" },
            { "mData": "browser", "sTitle": "Browser" },
            { "mData": "platform", "sTitle": "Platform" },
            { "mData": "version", "sTitle": "Version", "sClass": "center" },
            { "mData": "grade", "sTitle": "Grade", "sClass": "center" }
        ]
    });
  }
});
