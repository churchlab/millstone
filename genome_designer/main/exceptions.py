"""
Common Exception and Error classes.
"""


class InputError(Exception):
    """Exception raised for errors in the input.

    I assume this is what DBG wanted in d922023e8a1d4e7757bfa3e68aba3ace89f2a749.

    Attributes:
        msg  -- explanation of the error
    """
    def __init__(self, msg):
        self.msg = msg


class ValidationException(Exception):
    """Thrown when validation fails.
    """
    pass
