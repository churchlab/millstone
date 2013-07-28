/**
 * @fileoverview List of variants.
 */


gd.AlignmentCreateView = Backbone.View.extend({
  el: '#gd-page-container',


  /** Override. */
  initialize: function() {
    this.render();
  },


  /** Override. */
  render: function() {
    $('#gd-sidenav-link-alignments').addClass('active');

    this.refGenomeDataTable = new gd.DataTableComponent({
        el: $('#gd-datatable-ref_genome-hook'),
        objList: REF_GENOME_LIST_DATA['obj_list'],
        fieldConfig: REF_GENOME_LIST_DATA['field_config']
    });

    this.samplesDatatable = new gd.DataTableComponent({
        el: $('#gd-datatable-samples-hook'),
        objList: SAMPLES_LIST_DATA['obj_list'],
        fieldConfig: SAMPLES_LIST_DATA['field_config']
    });
  }
});
