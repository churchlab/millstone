/**
 * @fileoverview Simple utility functions.
 */


/**
 * Returns a new string instance that is this instance with the first letter
 * capitalized.
 */
String.prototype.capitalize = function() {
  return this.charAt(0).toUpperCase() + this.slice(1);
}
