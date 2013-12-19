"""
Object that stores the results of a filter being applied to a
Variant query and provides convenience methods for combining results.
"""

def metadata_default_dict_factory_fn():
    """Factory function that can be passed to defaultdict to provide default
    values for
    """
    return {'passing_sample_ids': set()}


class FilterEvalResult(object):
    """Wraps the result of evaluating a filter condition.

    Provides utility methods for combining results.
    """

    def __init__(self, variant_set, variant_id_to_metadata_dict={}):
        if not isinstance(variant_set, set):
            variant_set = set(variant_set)
        self.variant_set = variant_set
        self.variant_id_to_metadata_dict = variant_id_to_metadata_dict

    def __or__(self, other):
        return self.combine(other, '|')

    def __and__(self, other):
        return self.combine(other, '&')

    def combine(self, other, op_string):
        """Method that returns a new FilterEvalResult that is the combination
        of this one and other.

        Args:
            other: The FilterEvalResult to combine with.
            op_string: Either '&' or '|'.

        Returns:
            A new FilterEvalResult object.
        """
        assert isinstance(other, FilterEvalResult)
        assert op_string in ['&', '|']

        # Combine the Variant sets.
        if op_string == '&':
            new_variant_set = self.variant_set & other.variant_set
        elif op_string == '|':
            new_variant_set = self.variant_set | other.variant_set
        else:
            raise AssertionError("Unsupported op: %s" % op_string)

        # Build up the new metadata map.
        new_variant_id_to_metadata_dict = {}
        for variant in new_variant_set:
            merged_filter_metadata = {}

            self_filter_metadata = self.variant_id_to_metadata_dict.get(
                    variant.id, {})
            other_filter_metadata = other.variant_id_to_metadata_dict.get(
                    variant.id, {})

            # Merge passing sample ids.
            self_passing_genomes = self_filter_metadata.get(
                    'passing_sample_ids', set())
            other_passing_genomes = other_filter_metadata.get(
                    'passing_sample_ids', set())
            if op_string == '&':
                merged_filter_metadata['passing_sample_ids'] = (
                        self_passing_genomes & other_passing_genomes)
            else:
                merged_filter_metadata['passing_sample_ids'] = (
                        self_passing_genomes | other_passing_genomes)

            # Save the filter metadata.
            new_variant_id_to_metadata_dict[variant.id] = merged_filter_metadata

        return FilterEvalResult(new_variant_set,
                new_variant_id_to_metadata_dict)
