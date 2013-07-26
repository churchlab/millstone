#!/usr/bin/env python

"""
Setup script.
"""

import os

TOOLS_DIR = 'tools'

def setup():
    # Create tools dir if it doesn't exist.
    if not os.path.exists(TOOLS_DIR):
        os.mkdir(TOOLS_DIR)

    # Download the required tools from our Dropbox.

if __name__ == '__main__':
    setup()
