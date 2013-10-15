"""
Filter scope object.
"""

from main.models import ExperimentSample

FILTER_SCOPE__ALL = 'ALL'
FILTER_SCOPE__ANY = 'ANY'
FILTER_SCOPE__ONLY = 'ONLY'
VALID_FILTER_SCOPES = set([
    FILTER_SCOPE__ALL,
    FILTER_SCOPE__ANY,
    FILTER_SCOPE__ONLY
])

class FilterScope(object):
    """Represents the scope that a filter should be applied over.
    """

    def __init__(self, scope_type, sample_ids):
        """
        Args:
            sample_ids: Set of sample ids.
            scope_type: A scope in VALID_FILTER_SCOPES.
        """
        assert scope_type in VALID_FILTER_SCOPES, "Invalid scope type."

        self.sample_id_set = set(sample_ids)
        self.scope_type = scope_type


    def do_passing_samples_satisfy_scope(self, samples_passing_for_variant):
        """Returns a Boolean indicating whether the samples satisfy the scope.
        """
        if self.scope_type == FILTER_SCOPE__ALL:
            # All passing sample ids must be in the
            # scope set.
            intersection = (samples_passing_for_variant &
                    self.sample_id_set)
            return (intersection == self.sample_id_set)
        elif self.scope_type == FILTER_SCOPE__ANY:
            # At least one passing sample id must be in
            # the scope set.
            return len(samples_passing_for_variant &
                    self.sample_id_set) > 0
        elif self.scope_type == FILTER_SCOPE__ONLY:
            # The passing sample id set must be exactly
            # the scope set.
            return (samples_passing_for_variant ==
                    self.sample_id_set)

    @classmethod
    def parse_sample_ids(clazz, sample_id_string):
        """Turns a comma-separated list of sample uids or names into ids.
        """
        sample_uids_or_names = sample_id_string.split(',')
        sample_uids_or_names = [s.strip() for s in sample_uids_or_names]
        sample_ids = [ExperimentSample.objects.get(uid=uid).id for uid
                in sample_uids_or_names]
        return sample_ids
