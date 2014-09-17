#!/usr/bin/env python
"""
Basic tests for the fragmentation model
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import fragmentation_model as fm
import gene_cluster_library as gcl

def run_fragmentation_test ():
	"""Test and plot some fragmentation profiles.
	"""
	# We need to do some plotting
	import matplotlib.pyplot as plt
	# Some default values for the fragmentation data
	frag_mean = 250
	frag_sd = 75
	# Some example mRNA lengths
	mrna1_len = 500
	mrna2_len = 1000
	mrna3_len = 2000
	# Generate the profiles
	p1 = fm.frag_factor_profile(mrna1_len, frag_mean=frag_mean, frag_sd=frag_sd, no_end=False)
	p2 = fm.frag_factor_profile(mrna2_len, frag_mean=frag_mean, frag_sd=frag_sd, no_end=False)
	p3 = fm.frag_factor_profile(mrna3_len, frag_mean=frag_mean, frag_sd=frag_sd, no_end=True)
	# Plot the profiles
	fig = plt.figure(figsize=(8,3))
	ax = fig.add_subplot(1,1,1)
	ax.plot(range(mrna1_len), p1, 'r')
	ax.plot(range(mrna2_len), p2, 'g')
	ax.plot(range(mrna3_len), p3, 'b')
	ax.set_xlim([0,mrna3_len])
	ax.set_ylim([0,1.1])
	plt.tight_layout()
	plt.savefig('fragmentation_profile.pdf')

def run_fragmentation_variant_test ():
	# Load the Stata nif library data
	nifs = gcl.GeneClusterLibrary()
	nifs.load('../../data_sets/clean_nif_stata_library.txt')
	# Generate profiles
	fwd_profiles, rev_profiles = fm.variant_frag_factor_profile(nifs, '75', frag_mean=280, 
		                                                        frag_sd=70, bp_cutoff=10*70)
	# There should be 4 FWD and 3 REV profiles
	print 'Forward profiles = ', len(fwd_profiles) , '(Should equal 4)'
	print 'Reverse profiles = ', len(rev_profiles) , '(Should equal 3)'

# Test the functions
run_fragmentation_test()
run_fragmentation_variant_test()
