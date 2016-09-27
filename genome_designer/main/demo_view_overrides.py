"""View overrides for demo mode.
"""

from django.conf import settings
from django.contrib.auth import authenticate
from django.contrib.auth import login
from django.http import HttpResponseRedirect


def login_demo_account(request):
    new_user = authenticate(
            username=settings.DEMO_MODE__USERNAME,
            password=settings.DEMO_MODE__PASSWORD)
    login(request, new_user)
    redirect_url = request.GET.get('next', '/')
    return HttpResponseRedirect(redirect_url)
