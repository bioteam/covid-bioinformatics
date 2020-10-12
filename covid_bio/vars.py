#!/usr/bin/env python3

import os

'''
Create a default destination directory for all downloaded and 
created files, specify email address for Entrez.
'''

COV_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'COV')

EMAIL = 'self@organization.net'