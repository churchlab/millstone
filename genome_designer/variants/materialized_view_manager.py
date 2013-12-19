"""
Manages the Materialized view of the Variant data for filtering.
"""

from django.db import connection
from django.db import transaction

class AbstractMaterializedViewManager(object):
    """Base class for object acting as wrapper for a Postgresql materialized
    view.
    """

    def __init(self):
        self.view_table_name = None


    def create(self):
        """Creates the materialized view in the Postgresql DB.

        Child classes should implement.
        """
        # Drop the existing table.
        self.drop()

        # Delegate to child class.
        self.create_internal()


    def create_internal(self):
        """Creates the materialized view in the Postgresql DB.

        Child classes should implement.
        """
        pass


    def refresh(self):
        """Refreshes the view.
        """
        refresh_statement = 'REFRESH MATERIALIZED VIEW %s ' % (
                self.view_table_name)
        self.cursor.execute(refresh_statement)


    def drop(self):
        """Drops the materialized view in the Postgresql DB.
        """
        assert self.view_table_name
        drop_sql_statement = "DROP MATERIALIZED VIEW IF EXISTS %s" % (
                self.view_table_name,)
        self.cursor.execute(drop_sql_statement)
        transaction.commit_unless_managed()


SELECT_FIELDS = [
    'main_variant.id',
    'main_variant.position',
    'main_variant.chromosome',
    'main_experimentsample.id AS experiment_sample_id',
    'main_experimentsample.uid AS experiment_sample_uid'
]

class MeltedVariantMaterializedViewManager(AbstractMaterializedViewManager):
    """Interface for objects providing a wrapper for a Postgresql materialized
    view.
    """

    def __init__(self, reference_genome):
        self.reference_genome = reference_genome
        self.view_table_name = 'materialized_melted_variant'
        self.cursor = connection.cursor()

    def create_internal(self):
        """Override.
        """
        create_sql_statement = (
            'CREATE MATERIALIZED VIEW %s AS '
                'SELECT %s FROM main_variant '
                    'INNER JOIN main_variantcallercommondata ON (main_variant.id = main_variantcallercommondata.variant_id) '
                    'INNER JOIN main_variantevidence ON (main_variantcallercommondata.id = main_variantevidence.variant_caller_common_data_id) '
                    'INNER JOIN main_experimentsample ON (main_variantevidence.experiment_sample_id = main_experimentsample.id) '
                'WHERE (main_variant.reference_genome_id = %d) '
                'ORDER BY main_variant.position'
            % (self.view_table_name, ', '.join(SELECT_FIELDS), self.reference_genome.id)
        )
        self.cursor.execute(create_sql_statement)
        transaction.commit_unless_managed()
