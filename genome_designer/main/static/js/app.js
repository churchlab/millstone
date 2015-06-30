/**
 * The main app component.
 */


/**
 * For now, simply a router that reads the template var VIEW_TAG
 * and decides which view component to instantiate.
 */
gd.App = function() {

};


/** Runs the app. */
gd.App.prototype.run = function() {
  if (typeof VIEW_TAG === 'undefined') {
    return;
  }

  // Route depending on the VIEW_TAG in the template.
  switch(VIEW_TAG) {
    case 'PROJECT':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.ProjectView({'model': model});
      break;

    case 'REF_GENOME':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.RefGenomeView({'model': model});
      break;

    case 'REF_GENOME_LIST':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.RefGenomeListView({'model': model});
      break;

    case 'SAMPLE':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.SampleListView({'model': model});
      break;

    case 'ALIGNMENT':
      var model = new Backbone.Model(INIT_JS_DATA);
      var view = new gd.AlignmentView({'model': model});
      break;

    case 'ALIGNMENT_LIST':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.AlignmentListView({'model': model});
      break;

    case 'ALIGNMENT_CREATE':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.AlignmentCreateView({'model': model});
      break;

    case 'VARIANT_SET':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.VariantSetView({'model': model});
      break;

    case 'VARIANT_SET_LIST':
      var model = new Backbone.Model(INIT_JS_DATA.entity);
      var view = new gd.VariantSetListView({'model': model});
      break;

    case 'GENE':
      var view = new gd.GeneView();
      break;

    case 'GOTERM':
      var view = new gd.GotermView();
      break;

    case 'ANALYZE':
      var model = new Backbone.Model(INIT_JS_DATA);
      var view = new gd.TabAnalyzeBaseView({'model': model});
      break;
  }
};


$(document).ready(function() {
  var app = new gd.App();
  app.run();
  Backbone.history.start({pushState: true});
});
