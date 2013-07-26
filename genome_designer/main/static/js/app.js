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
    case 'REF_GENOME':
      var view = new gd.RefGenomeView();
      break;

    case 'SAMPLE':
      var view = new gd.SampleView();
      break;

    case 'ALIGNMENT':
      var view = new gd.AlignmentView();
      break;
    
    case 'ALIGNMENT_CREATE':
      var view = new gd.AlignmentCreateView();
      break;
    

    case 'VARIANT_SET':
      var view = new gd.VariantSetView();
      break;

    case 'VARIANT':
      var view = new gd.VariantView();
      break;

    case 'GENE':
      var view = new gd.GeneView();
      break;

    case 'GOTERM':
      var view = new gd.GotermView();
      break;
  }
};


$(document).ready(function() {
  var app = new gd.App();
  app.run();
});
