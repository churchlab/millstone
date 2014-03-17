"""
Utility objects and functions for debugging.
"""


class FakeException(Exception):
    """Placeholder exception for swapping into the except clause of a
    try-except clause so that errors are thrown.

    For example, this is useful in xhr_handlers.get_variant_list() where I have
    a try-except that returns a 404 to the user, rather than a 500. However,
    for debugging, we often want to see the full stacktrace.
    """
    pass
