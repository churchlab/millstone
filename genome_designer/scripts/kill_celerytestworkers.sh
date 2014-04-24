#!/bin/sh
# Script to kill celery test workers.
kill $(ps aux | grep '[p]ython.*celerytestworker' | awk '{print $2}')
