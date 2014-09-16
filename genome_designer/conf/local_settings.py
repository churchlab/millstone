# Add settings overrides here.
# Changes to this file are ignored by git.

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'gdv2db',
        'USER': 'gdv2dbuser',
        'PASSWORD': 'g3n3d3z',
        'HOST': 'localhost',
        'PORT': '',
        'OS_USER': 'dbgoodman'
    }

    # Uncomment to debug with Sqlite.
    # 'default': {
    #     'ENGINE': 'django.db.backends.sqlite3',
    #     'NAME': 'temp.db',
    #     'USER': '',
    #     'PASSWORD': '',
    #     'HOST': '',
    #     'PORT': '',
    # }
}