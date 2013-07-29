/**
 * Decorates the toggle-able part of the nav component when present.
 */


$(document).ready(function() {
  $('.tree-toggle').click(function () {
    $(this).parent().children('ul.tree').toggle(200);
  });
});
