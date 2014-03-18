gd.Spinner = function() {
  var opts = {
      lines: 13, // The number of lines to draw
      length: 20, // The length of each line
      width: 10, // The line thickness
      radius: 30, // The radius of the inner circle
      corners: 1, // Corner roundness (0..1)
      rotate: 0, // The rotation offset
      direction: 1, // 1: clockwise, -1: counterclockwise
      color: '#000', // #rgb or #rrggbb or array of colors
      speed: 1, // Rounds per second
      trail: 60, // Afterglow percentage
      shadow: false, // Whether to render a shadow
      hwaccel: false, // Whether to use hardware acceleration
      className: 'spinner', // The CSS class to assign to the spinner
      zIndex: 2e9, // The z-index (defaults to 2000000000)
      top: ($(window).height() / 2) + 'px', // Top position relative to parent in px
      left: ($(window).width() / 2) + 'px'// Left position relative to parent in px

  };
  this.spinnerDelegate = new Spinner(opts);
};


/**
 * Attaches the spinner to the target, or body by default, and starts
 * spinning.
 */
gd.Spinner.prototype.spin = function(opt_target) {
  var target = opt_target || document.body;
  console.log(this.spinnerDelegate);
  this.spinnerDelegate.spin(target);
};


/** Stops the spinner. */
gd.Spinner.prototype.stop = function() {
  this.spinnerDelegate.stop();
};
