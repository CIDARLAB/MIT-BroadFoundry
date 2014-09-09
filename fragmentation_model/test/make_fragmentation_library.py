#!/usr/bin/env python
"""
Basic tests for the fragmentation model
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import fragmentation_model as fm

# Create and save the library
fm.save_frag_profile_library('test_frag_library.txt', 1, 200, frag_mean=280, frag_sd=70, bp_cutoff=None)

# Load the library
profiles_fixed, profiles_open = fm.load_frag_profile_library ('test_frag_library.txt')
