/**
 * @fileoverview Base object for views related to analysis.
 *
 * TODO: Refactor to use gd.VariantsTableComponent.
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
    // Decorate sidebar.
    $('.gd-sidebar > div.list-group > select').selectpicker()

    this.registerManualListeners();

    if (this.model.has('alignmentGroup')) {
      this.handleAlignmentGroupSelect(undefined,
          this.model.get('alignmentGroup').uid,
          this.model.get('refGenome').uid);
    } else if (this.model.has('refGenome')) {
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
    $('#gd-analyze-select-ref-genome').selectpicker();

    $('#gd-analyze-select-ag').change(
        _.bind(this.handleAlignmentGroupSelect, this));
    $('#gd-analyze-select-ag').selectpicker('hide');

    $('#gd-analyze-select-search-entity').change(
        _.bind(this.handleSearchEntitySelect, this));
    $('#gd-analyze-select-search-entity').selectpicker('hide');
  },


  /**
   * Handles a Reference Genome being selected from the dropdown menu.
   *
   * This method reveals the AlignmentGroup dropdown.
   *
   * @param {Event} e Event object if this is being called in response to an
   *     event. Ignored if opt_refGenomeUid is provided.
   * @param {string=} opt_refGenomeUid Optional ref genome string explicitly
   *     provided.
   */
  handleRefGenomeSelect: function(e, opt_refGenomeUid) {
    var refGenomeUid = opt_refGenomeUid || $(e.target).val();
    this.model.set('refGenomeUid', refGenomeUid);

    // Destroy the current subview.
    this.destroySubview();

    // Initialize the AlignmentGroup dropdown.
    $('#gd-analyze-select-ag').empty();
    $('#gd-analyze-select-ag').append(
        '<option selected disabled>Alignment Group</option>');

    // Hide the entity select if showing.
    //$('#gd-analyze-select-search-entity').hide();

    // Populate the AlignmentGroup select with the relevant options and fade
    // in.
    var requestData = {
      'refGenomeUid': this.model.get('refGenomeUid')
    };

    $.get('/_/alignmentgroups', requestData,
        _.bind(function(response) {

          // ReferenceGenome has no AlignmentGroups.
          if (!response.length > 0) {
            $('#gd-analyze-select-ag').append(
                '<option selected disabled>No Alignments</option>');
            $('#gd-analyze-select-ag').selectpicker('refresh');
            $('#gd-analyze-select-search-entity').selectpicker('hide');
            return;
          }

          // Otherwise, add option for each AlignmentGroup.
          _.each(response, function(alignmentGroup) {
            $('#gd-analyze-select-ag').append(
                '<option value=' + alignmentGroup.uid + '>' +
                  alignmentGroup.label +
                '</option>');
          })

          $('#gd-analyze-select-ag').selectpicker('refresh');

        }, this));

    $('#gd-analyze-select-ag').selectpicker('show');
  },


  /**
   * Handles an AlignmentGroup being selected from the dropdown menu.
   *
   * This method does an initial fetch of the data, and (re)-draws the
   * DataTable component. Subsequent requests given the same filter are
   * handled through a pagination dance between the DataTables component
   * and server-side support.
   *
   * @param {Event} e Event object if this is being called in response to an
   *     event. Ignored if opt_refGenomeUid is provided.
   * @param {string=} opt_alignmentGroupUid Optional AlignmentGroup uid
   *     explicitly provided.
   * @param {string=} opt_refGenomeUid Optional ReferenceGenom uid
   *     explicitly provided.
   */
  handleAlignmentGroupSelect: function(e, opt_alignmentGroupUid,
      opt_refGenomeUid) {
    var alignmentGroupUid = opt_alignmentGroupUid || $(e.target).val();
    this.model.set('alignmentGroupUid', alignmentGroupUid);

    // Manually setting this if not already set during previous select.
    if (opt_refGenomeUid) {
      this.model.set('refGenomeUid', opt_refGenomeUid);
    }

    // Display the AG dropdown if not already shown
    $('#gd-analyze-select-ag').selectpicker('show');

    // Determine the subview.
    var subViewKey = 'variants'; // Default
    if (this.model.has('subView')) {
      subViewKey = this.model.get('subView');
    }

    // Update the url.
    this.router.navOnAlignmentGroupSelect(alignmentGroupUid, subViewKey);

    // Show the entity select (i.e. Variants, VariantSets, etc.)
    $('#gd-analyze-select-search-entity').selectpicker('show')

    // Determine the sub-view to show from the key.
    var subView = gd.TabAnalyzeBaseView.SUBVIEW_TYPE_TO_VIEW_CLASS[
        subViewKey.toUpperCase()];
    this.updateSubview(subView);
  },


  /**
   * Handles the entity select being triggered.
   *
   * @param {Event} e Event object if this is being called in response to an
   *     event.
   */
  handleSearchEntitySelect: function(e) {
    // Get the view key from DOM.
    var subViewOptionValue = $(e.target).val();

    // Update the url.
    var subViewUrlToken = subViewOptionValue.toLowerCase();
    var alignmentGroupUid = this.model.get('alignmentGroupUid');
    this.router.navOnAlignmentGroupSelect(alignmentGroupUid, subViewUrlToken);

    // Update the model subview key.
    this.model.set('subView' ,subViewUrlToken);

    // Update the view.
    var newSubViewType = gd.TabAnalyzeBaseView.SUBVIEW_TYPE_TO_VIEW_CLASS[
        subViewOptionValue];
    this.updateSubview(newSubViewType);
  },


  /** Draws or replaces the current subview with the new selected subview. */
  updateSubview: function(newSubviewType) {
    this.destroySubview();

    this.currentSubView = new newSubviewType({model: this.model});
    this.listenTo(this.currentSubView, 'NAVIGATE',
        _.bind(this.handleSubviewNav, this));
  },


  /** Tear out the current subview. */
  destroySubview: function() {
    // Get rid of the current view.
    $('.gd-table-controls').empty();
    if (this.currentSubView) {
      this.currentSubView.destroy();
    }
    this.currentSubView = null;

    // Remove anything left there.
    $('#gd-datatable-hook').empty();
  },


  /** Handles a nav event triggered by a subview. */
  handleSubviewNav: function(navParams) {
    this.router.navigate(navParams.path);
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
