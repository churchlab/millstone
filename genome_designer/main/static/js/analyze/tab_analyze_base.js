/**
 * @fileoverview Base object for views related to analysis.
 */

gd.TabAnalyzeBaseView = Backbone.View.extend(
/** Prototype properties */
{
  /** Element that Backbone attches this component to. */
  el: '#gd-page-container',


  /** Backbone initialize method. */
  initialize: function() {
    // Router object for synchronizing history.
    this.router = new gd.TabAnalyzeRouter({
        projectUid: this.model.get('project').uid
    });

    // Identify the current subview that is visible.
    this.currentSubView = null;

    this.render();
  },


  /** Override. */
  render: function() {
    this.registerManualListeners();

    if (this.model.has('refGenome')) {
      this.handleRefGenomeSelect(undefined, this.model.get('refGenome').uid);
    }
  },


  /**
   * Registers listeners that can't be registered using the Backbone syntax
   * sugar 'events' property.
   */
  registerManualListeners: function() {
    $('#gd-analyze-select-ref-genome').change(
        _.bind(this.handleRefGenomeSelect, this));
    $('#gd-analyze-select-search-entity').change(
        _.bind(this.handleSearchEntitySelect, this));
  },


  /**
   * Handles a Reference Genome being selected from the dropdown menu.
   *
   * This method does an initial fetch of the data, and (re)-draws the
   * DataTable component. Subsequent requests given the same filter are
   * handled through a pagination dance between the DataTables component
   * and server-side support.
   *
   * @param {Event} e Event object if this is being called in response to an
   *     event. Ignored if opt_refGenomeUid is provided.
   * @param {string=} opt_refGenomeUid Optional ref genome string explicitly
   *     provided.
   */
  handleRefGenomeSelect: function(e, opt_refGenomeUid) {
    var refGenomeUid = opt_refGenomeUid || $(e.target).val();
    this.model.set('refGenomeUid', refGenomeUid);

    // Update the url.
    this.router.navOnRefGenomeSelect(refGenomeUid);

    // Show the entity select (i.e. Variants, VariantSets, etc.)
    // TODO: Implement support for switching among entities.
    $('#gd-analyze-select-search-entity').fadeIn();

    this.updateSubview(
        gd.TabAnalyzeBaseView.SUBVIEW_TYPE_TO_VIEW_CLASS.VARIANTS);
  },


  /**
   * Handles the entity select being triggered.
   *
   * @param {Event} e Event object if this is being called in response to an
   *     event.
   */
  handleSearchEntitySelect: function(e) {
    var newSubviewType = gd.TabAnalyzeBaseView.SUBVIEW_TYPE_TO_VIEW_CLASS[
        $(e.target).val()];
    this.updateSubview(newSubviewType);
  },


  /** Draws or replaces the current subview with the new selected subview. */
  updateSubview: function(newSubviewType) {
    // Get rid of the current view.
    $('#gd-analyze-subview-controls-hook').empty();
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }
    this.datatableComponent = null;

    // Remove anything left there.
    $('#gd-datatable-hook').empty();

    this.currentSubView = new newSubviewType({model: this.model});
  }
},

/** Static properties */
{
  /** Mape from subview to the class that decorates that view. */
  SUBVIEW_TYPE_TO_VIEW_CLASS: {
      GENES: gd.TabAnalyzeSubviewGenes,
      SETS: gd.TabAnalyzeSubviewSets,
      VARIANTS: gd.TabAnalyzeSubviewVariants
  }
}
);
