/**
 * @fileoverview Simple utility functions.
 */

// Namespace.
gd.Util = {};


// Enable Select-Picker.
$('.selectpicker').selectpicker();

/**
 * Returns a new string instance that is this instance with the first letter
 * capitalized.
 */
String.prototype.capitalize = function() {
  return this.charAt(0).toUpperCase() + this.slice(1);
};


/**
 * Modify the URL query params at the given key.
 *
 * Copied and modified from: http://stackoverflow.com/a/6021027/397516
 * No comment on bugs at this point.
 *
 * @return {string} Updated uri with the query string modified.
 */
gd.Util.updateQueryStringParameter = function(uri, key, value) {
  var re = new RegExp("([?|&])" + key + "=.*?(&|$)", "i");
  var separator = uri.indexOf('?') !== -1 ? "&" : "?";
  if (uri.match(re)) {
    return uri.replace(re, '$1' + key + "=" + value + '$2');
  }
  else {
    return uri + separator + key + "=" + value;
  }
};



/** Returns the current full path, including query params. */
gd.Util.getFullPath = function() {
  // Use the jQuery location wrapper object.
  return window.location.pathname + window.location.search;
};


/**
 * Returns an object representing the current query string.
 *
 * NOTE: This seems like something that should be in some javascript library
 * but I can't find one, hence I wrote this potentially buggy implementation.
 */
gd.Util.getQueryStringAsObject = function() {
  var paramObject = {};

  var maybeQueryString = window.location.search;
  if (maybeQueryString) {
    maybeQueryString = maybeQueryString.slice(1, maybeQueryString.length);
    var paramList = maybeQueryString.split('&');
    _.each(paramList, function(param) {
        // Find location of first equal sign and split at that point.
        var dividerIndex = param.indexOf('=');
        var key = param.slice(0, dividerIndex);
        var value = param.slice(dividerIndex + 1, param.length);
        paramObject[decodeURIComponent(key)] = decodeURIComponent(value);
    });
  }

  return paramObject;
};
