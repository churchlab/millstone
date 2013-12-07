#!/usr/bin/env python

"""
Convenience script for parsing profiler output.
"""

import hotshot.stats
import sys

# Read the file
stats = hotshot.stats.load(sys.argv[1])

# Uncomment to strip directories.
#stats.strip_dirs()

# This sort should reveal time spent per function.
stats.sort_stats('cumulative', 'calls')

# Write top 40 stats
stats.print_stats(40)
