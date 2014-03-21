"""
Custom Django middleware.

Middleware is Django's interface to allow actions to wrap requests
and views.

See: https://docs.djangoproject.com/en/dev/topics/http/middleware/
"""

from django.conf import settings
from django.http import HttpResponseBadRequest

from conf import demo_settings


class DisabledInDemoModeMiddleware(object):
    """Disable all views except those that are explicitly allowed.
    """

    def __init__(self):
        """Constructor.
        """
        # Build up a list of safe views.
        self.safe_views = []
        if hasattr(demo_settings, 'DEMO_SAFE_VIEWS'):
            for view_path in demo_settings.DEMO_SAFE_VIEWS:
                view = self._get_view(view_path)
                self.safe_views.append(view)
        print self.safe_views

    def _get_view(self, view_path):
        """Get the View from the path to help build the list of safe views.
        """
        right_most_dot = view_path.rfind('.')
        module_path, view_name = (
                view_path[:right_most_dot], view_path[right_most_dot + 1:])
        module = __import__(module_path, globals(), locals(), [view_name])
        return getattr(module, view_name)

    def process_view(self, request, view_func, *view_args, **view_kwargs):
        """Override to implement Middleware functionality.
        """
        # Nothing to do when not demo mode.
        if not settings.DEMO_MODE:
            return None

        if view_func in self.safe_views:
            return None # continue handling request
        return HttpResponseBadRequest() # interrupt
