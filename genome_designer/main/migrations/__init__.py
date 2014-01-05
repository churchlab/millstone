# This code prevents south logging from showing up when a test fails.
import logging
south_logger=logging.getLogger('south')
south_logger.setLevel(logging.INFO)
