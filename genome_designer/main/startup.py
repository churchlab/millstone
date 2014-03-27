"""Actions to run at server startup.
"""

from django.db import connection
from django.db import transaction


def run():
    """Call this from manage.py or tests.
    """
    _add_custom_mult_agg_function()


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
