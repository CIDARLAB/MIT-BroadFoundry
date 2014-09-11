#!/usr/bin/env python
"""Functions to correct RNA-seq based traces using known fragmentation patterns.
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

# Load the Stata nif library data (ignore ones with sequencing errors => "clean")
nifs = gcl.GeneClusterLibrary()
nifs.load('../../data_sets/clean_nif_stata_library.txt')

# Generate expected profiles
fwd_profiles, rev_profiles = fm.variant_frag_factor_profile(nifs, '75', frag_mean=280, 
	                                                        frag_sd=70, bp_cutoff=10*70)

# For each TU attempt to fit to real trace data




