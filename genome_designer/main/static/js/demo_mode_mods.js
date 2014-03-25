/**
 * @fileoverview Modifications to make demo line work better.
 *
 * For example, we want to gracefully handle failures due to server endpoints
 * that are restricted in demo mode.
 */


// Register a listener for 403 errors, which only failed conneciton to
// demo pages should throw.
// TODO: Confirm other cases don't throw this.
$(document).ajaxError(function(event, jqXHR, ajaxSettings, thrownError) {
  if (jqXHR.status == 403) {
    alert('Sorry, this feature is restricted in demo mode.');
  }
});
