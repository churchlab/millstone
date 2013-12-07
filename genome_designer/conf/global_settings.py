# Django settings for genome_designer project.

import os

from django.conf import global_settings

# EntrezGene wants an email to use it's API.
EMAIL = "millstone_user@gmail.com"

# The absolute path of the settings.py file's directory.
# Useful for settings that require absolute paths like templates.
PWD = os.path.join(os.path.dirname(os.path.realpath(__file__ )), '../')

# Absolute path to the third-party tools dir where setup.py stores
# downloaded tools that are used internally.
TOOLS_DIR = os.path.join(PWD, 'tools')

DEBUG = True

TEMPLATE_DEBUG = DEBUG

# Default True, requiring celery server to be running.
# Override this as False in local_settings.py to force synchronous behavior.
DEBUG_CONCURRENT = True

ADMINS = (
    # ('Your Name', 'your_email@example.com'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'temp.db',
        'USER': '',                      # Not used with sqlite3.
        'PASSWORD': '',                  # Not used with sqlite3.
        'HOST': '',                      # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '',                      # Set to empty string for default. Not used with sqlite3.
    }
}

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/{{ docs_version }}/ref/settings/#allowed-hosts
ALLOWED_HOSTS = []

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'America/New_York'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = 'temp_data'

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = ''

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = ''

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = '/static/'

# URL prefix for admin static files -- CSS, JavaScript and images.
# Make sure to use a trailing slash.
# Examples: "http://foo.com/static/admin/", "/static/admin/".
ADMIN_MEDIA_PREFIX = '/static/admin/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = '*2s+5aalcl*0relwgx(4h9tguj8xt@j$fpr_&utx1!m%j70e()'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
)

ROOT_URLCONF = 'genome_designer.urls'

TEMPLATE_DIRS = (
    os.path.join(PWD, 'main/templates')
)

INSTALLED_APPS = (
    # django built-ins
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',

    # async queue
    'djcelery',

    # third-party apps
    'registration',

    # our apps
    'main',

    # Testing
    'django_nose'
)


LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,

    'formatters': {
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
    },

    'handlers': {
        'console':{
            'level':'DEBUG',
            'class':'logging.StreamHandler',
            'formatter': 'simple'
        },
        'file': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': 'default.log', # override in local_settings.py
            'formatter': 'simple'
        },
    },

    'loggers': {
        'django': {
            'handlers':['console'],
            'level':'INFO',
            'propagate': True,
        },

        's3': {
            'handlers':['console'],
            'level':'INFO',
            'propagate': True,
        },

        # Uncomment to see SQL logs on the console.
        'django.db.backends': {
            'handlers':['file'],
            'level':'DEBUG',
            'propagate': False,
        },
    }
}


###############################################################################
# Custom Config
###############################################################################

# We use a separate UserProfile from the build-in Django User model so that we
# have the option of extending it.
AUTH_PROFILE_MODULE = 'main.UserProfile'

# Custom template context processors.
TEMPLATE_CONTEXT_PROCESSORS = global_settings.TEMPLATE_CONTEXT_PROCESSORS + (
    'main.context_processors.common_data',
    'main.context_processors.aws_settings'
)

###############################################################################
# Registration / Accounts
###############################################################################

# django-registration
# One-week activation window.
ACCOUNT_ACTIVATION_DAYS = 7

LOGIN_REDIRECT_URL = '/'


###############################################################################
# Django Celery - async task queue management
###############################################################################

import djcelery
djcelery.setup_loader()

# RabbitMQ settings
BROKER_URL = 'amqp://guest:guest@localhost:5672//'


###############################################################################
# External tools
###############################################################################

TOOLS_DIR = 'tools'
BASH_PATH = '/bin/bash'

###############################################################################
# JBrowse
###############################################################################

# Root of the JBrowse installation.
JBROWSE_ROOT = os.path.abspath(os.path.join(PWD, '../jbrowse'))

# The location of the JBrowse scripts (perl scripts, ugh)
JBROWSE_BIN_PATH = os.path.abspath(os.path.join(PWD, '../jbrowse/bin'))

# Name of the symlink from within JBrowse to the data dir.
# See JBROWSE_DATA_URL_ROOT description below.
JBROWSE_DATA_SYMLINK_NAME = 'gd_data'

# Full path to the JBrowse symlink (links back to the app data dir).
JBROWSE_DATA_SYMLINK_PATH = os.path.join(JBROWSE_ROOT,
        JBROWSE_DATA_SYMLINK_NAME)

# The url root to data that JBrowse displays.
# The app admin should create a symlink from the actual data root to this
# location inside of the jbrowse/ dir. For example, the way we display bam
# files is configuring the trackList.json file with a track with the following
# key-value: "urlTemplate" : "/jbrowse/gd_data/users/8fc1f831/projects/58a62c7d/genomes/8dc829ec/align.bam"
JBROWSE_DATA_URL_ROOT= '/jbrowse/' + JBROWSE_DATA_SYMLINK_NAME + '/'

# Set to True if you want to force JBrowse links to be from localhost and thus
# go through nginx. Default of False causes JBrowse to serve through Django in
# dev.
DEBUG_FORCE_JBROWSE_NGINX = False

###############################################################################
# Snp Calling
###############################################################################

# Path to snpeff java jar.
SNPEFF_JAR_PATH = os.path.abspath(os.path.join(PWD, 'tools','snpEff',
        'snpEff.jar'))
# Path to snpeff config template.
SNPEFF_CFG_TEMPLATE_PATH = os.path.join(PWD, 'main',
        'templates','snpeff.tmpl.config')

# Upstream/downstream interval where SNPs nearby a gene are tagged. Needs to
#   be smaller than default for bacterial genomes.
SNPEFF_UD_INTERVAL_LENGTH = 50

# SNPEff can be multithreaded; probably makes sense to set this to the same as
# a more general CPU thread value.
SNPEFF_THREADS = 2

# If we're debugging snpeff, print the output
SNPEFF_BUILD_DEBUG = True

# Names of SnpEff summary files, which we want to delete after running.
SNPEFF_SUMMARY_FILES = ['snpEff_genes.txt', 'snpEff_summary.html']

###############################################################################
# S3
###############################################################################

# Allows user to create an S3 backed project. 
S3_ENABLED = True

# Don't perform any API call that changes anything on S3. 
S3_DRY_RUN = False

# Get them from https://console.aws.amazon.com/iam/home?#security_credential
AWS_CLIENT_SECRET_KEY = '69kDz9o8VzD1Avf981R40Yf+zioz2dwnY5g94UK4'
AWS_SERVER_PUBLIC_KEY = 'AKIAJERCJ7K6W7F7VE7A'
AWS_SERVER_SECRET_KEY = '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'

# Name of S3 bucket to which all files will be uploaded
S3_BUCKET = 'genome-designer-upload'

# Run S3 tests on S3_TEST_BUCKET; fine to use S3_BUCKET for S3_TEST_BUCKET, for
# all test operations will run inside s3://S3_TEST_BUCKET/__tests__
S3_TEST_BUCKET = S3_BUCKET

# Maximum file size for user upload
S3_FILE_MAX_SIZE = 1024 ** 3  # 1GB

###############################################################################
# Testing
###############################################################################

TEST_RUNNER = 'test_suite_runner.CustomTestSuiteRunner'

TEST_FILESYSTEM_DIR = os.path.join(PWD, 'temp_test_data')

TEST_S3 = S3_ENABLED

###############################################################################
# Profiling
###############################################################################

# Directory where profiler logs will be stored. See README.md.
PROFILE_LOG_BASE = None
