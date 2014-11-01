# Add settings overrides here.
# Changes to this file are ignored by git.

# Uncomment the following if running Postgres under different username
# (i.e. for Postgres.app on OSX)
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'gdv2db',
        'USER': 'gdv2dbuser',
        'PASSWORD': 'g3n3d3z',
        'HOST': 'localhost',
        'PORT': '',
        'OS_USER': 'YOUR_USERNAME'
    }
}