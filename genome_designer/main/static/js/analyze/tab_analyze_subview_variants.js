/**
 * @fileoverview Sub-view displaying Variants.
 */


gd.TabAnalyzeSubviewVariants = gd.TabAnalyzeSubviewAbstractBase.extend(
{
  /** @override */
  initialize: function() {
    this.variantsTableComponent = null;

    // List of variant sets belonging to the user.
    this.variantSetList = null;

    // Map of variant keys, used to draw the fields dropdown.
    this.variantKeyMap = null;

    // Keys that should be visible.
    // When null, we are just showing the default.
    this.visibleKeyNames = null;

    // Check whether a filter is set on the view.
    var paramObject = gd.Util.getQueryStringAsObject();
    var maybeQueryString = window.location.search;
    if ('filter' in paramObject) {
      var filterString = paramObject['filter'];
      this.model.set('filterString', filterString);
    } else {
      this.model.set('filterString', '');
    }

    // Set the view to cast as default, unless melt=1 param in url.
    this.model.set('isMelted',
        'melt' in paramObject && paramObject['melt'] == '1');

    this.variantsTableComponent = new gd.VariantsTableComponent({
      model: this.model
    });

    this.setupListeners();
  },

  setupListeners: function() {
    // Propagate NAVIGATE event.
    // TODO: Is there more elegant way to do this?
    this.listenTo(this.variantsTableComponent, 'NAVIGATE', function (e) {
      this.trigger('NAVIGATE', e);
    });
  },

  /** @override */
  destroy: function() {
    if (this.variantsTableComponent) {
      this.variantsTableComponent.destroy()
    }
    this.variantsTableComponent = null;
  }
});
