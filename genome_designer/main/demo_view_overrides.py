"""View overrides for demo mode.
"""

from django.conf import settings

if settings.DEMO_MODE:
    from django.contrib.auth import authenticate
    from django.contrib.auth import login
    from django.http import HttpResponseRedirect

    def login_demo_account(request):
        new_user = authenticate(username='gmcdev',
                password='g3n3d3z')
        login(request, new_user)
        return HttpResponseRedirect("/")
