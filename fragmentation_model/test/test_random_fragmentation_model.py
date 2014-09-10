#!/usr/bin/env python
"""
Basic tests for the random fragmentation model to validate analytic solution
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import fragmentation_model as fm
import matplotlib.pyplot as plt

if False:
	profile = fm.frag_factor_profile_random_prob_model(700, 280, 70, mrna_count=10000, 
		                                               bp_frag_prob=1.0/280.0, verbose=True)
	fig = plt.figure(figsize=(8,3))
	ax = fig.add_subplot(1,1,1)
	ax.plot(range(700), profile, 'b')
	ax.set_xlim([0,700])
	ax.set_ylim([0,profile.max()*1.1])
	plt.tight_layout()
	plt.savefig('random_prob_fragmentation_profile.pdf')

if True:
	profile = fm.frag_factor_profile_random_break_model(1000, 280, 70, mrna_count=100000, 
		                                                max_frags=20, verbose=True)
	fig = plt.figure(figsize=(8,3))
	ax = fig.add_subplot(1,1,1)
	ax.plot(range(700), profile, 'b')
	ax.set_xlim([0,700])
	ax.set_ylim([0,profile.max()*1.1])
	plt.tight_layout()
	plt.savefig('random_break_fragmentation_profile.pdf')

