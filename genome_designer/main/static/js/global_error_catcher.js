/**
 * Catches any 500 errors and shows an alert to the user with a link to report
 * the bug. Otherwise, the error is hidden in the console and most users
 * won't know what's going on.
 */

 $(document).ajaxError(function (e, xhr, options) {
  var url = e.currentTarget.URL;
  if (xhr.status == 500) {
      alert(
          'Server Error. Please file an issue at ' +
          'https://github.com/churchlab/millstone/issues/new. ' +
          'Describe what you were trying to do, and include the following ' +
          'offending URL: ' + url);
  }
});
