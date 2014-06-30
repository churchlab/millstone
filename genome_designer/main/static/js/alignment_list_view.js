/**
 * @fileoverview List of alignments.
 */


gd.AlignmentListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
    this.refreshTable();
  },

  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.redrawDatatable();
  },

  refreshTable: function() {
    this.refreshInterval = window.setInterval(_.bind(this.redrawDatatable,this), 5000);
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-alignment_list_view-datatable-hook'),
        serverTarget: '/_/alignmentgroups',
        controlsTemplate: '/_/templates/alignment_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  decorateControls: function() {
    var dict_list = this.datatableComponent.getDataDictList();
    var stopRefreshing = true;

    for(var i = 0; i < dict_list.length; i++) {

      if(dict_list[i]['run_time']!='Not running') {
        stopRefreshing = false;
        break;
      }

    }

    if(stopRefreshing) {
      window.clearInterval(this.refreshInterval);
    }

    this.controlsComponent = new gd.AlignmentListControlsComponent({
      el: '#gd-sample-list-view-datatable-hook-control',
      datatableComponent: this.datatableComponent
    });

    this.listenTo(this.controlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  }
});
