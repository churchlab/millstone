from conf.global_settings import *
from conf.local_settings import *

DATABASES['default'] = {
    'ENGINE': 'django.db.backends.postgresql_psycopg2',
    'NAME': 'gdv2db',
    'USER': 'gdv2dbuser',
    'PASSWORD': 'g3n3d3z',
    'HOST': 'localhost',
    'PORT': '',
}
