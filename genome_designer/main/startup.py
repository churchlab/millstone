"""Actions to run at server startup.
"""

import subprocess
import sys

from django.db import connection
from django.db import transaction
from south.migration import Migrations
from south.models import MigrationHistory


def run():
    """Call this from manage.py or tests.
    """
    _add_custom_mult_agg_function()

    _check_environment()

    # TODO: This breaks test_pipeline.py. Why?
    # _check_migrations_applied()


def _add_custom_mult_agg_function():
    """Make sure the Postgresql database has a custom function array_agg_mult.

    NOTE: Figured out the raw sql query by running psql with -E flag
    and then calling \df to list functions. The -E flag causes the internal
    raw sql of the commands to be shown.
    """
    cursor = connection.cursor()

    cursor.execute(
            'SELECT p.proname '
            'FROM pg_catalog.pg_proc p '
            'WHERE p.proname=\'array_agg_mult\''
    )
    mult_agg_exists = bool(cursor.fetchone())
    if not mult_agg_exists:
        cursor.execute(
                'CREATE AGGREGATE array_agg_mult (anyarray)  ('
                '    SFUNC     = array_cat'
                '   ,STYPE     = anyarray'
                '   ,INITCOND  = \'{}\''
                ');'
        )
        transaction.commit_unless_managed()


def _check_migrations_applied():
    """Checks that all south migrations have been applied.
    """
    APP_NAME = 'main'
    all_migrations = Migrations(APP_NAME)
    applied_migrations = [migration.get_migration() for migration in
            MigrationHistory.objects.filter(app_name=APP_NAME)]
    not_applied = set(all_migrations) - set(applied_migrations)
    if len(not_applied):
        raise AssertionError(
                "Database migration required. "
                "Please run `./manage.py migrate main`.\n"
                "Applied: {applied}\n"
                "Missing: {not_applied}\n".format(
                        applied=applied_migrations,
                        not_applied=not_applied))


def _check_environment():
    """Checks for software and tools in environment.

    NOTE: In progress.
    """
    # Check for java.
    try:
        subprocess.check_output(
                ["java", "-version"], stderr=subprocess.STDOUT)
    except OSError:
        raise AssertionError("Startup Error: java not found in environment.")
