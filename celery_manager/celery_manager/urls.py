from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
import os.path

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'^$', 'celery_manager.views.home', name='home'),
    url(r'^test$', 'celery_manager.views.test', name='test'),
    url(r'^save$', 'celery_manager.views.save', name='save'),
    url(r'^status$', 'celery_manager.views.status', name='status'),
) + static(settings.STATIC_URL, document_root=os.path.join(settings.BASE_DIR, "static"))
