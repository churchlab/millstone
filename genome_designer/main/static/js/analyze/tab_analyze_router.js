/**
 * @fileoverview Object that syncs the URL with the state of the Analyze tab.
 */


gd.TabAnalyzeRouter = Backbone.Router.extend({


  /** Override. */
  initialize: function(options) {
    this.projectUid = options.projectUid;

    this.rootUrl = '/projects/' + this.projectUid + '/analyze';
  },


  /** Navigates to the ref genome uid. */
  navOnAlignmentGroupSelect: function(alignmentGroupUid, subViewKey) {
    gd.Spinner.globalClear();
    this.navigate(this.rootUrl + '/' + alignmentGroupUid + '/' + subViewKey);
  }
});
