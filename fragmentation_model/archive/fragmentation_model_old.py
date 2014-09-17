#!/usr/bin/env python
"""
Fragmentation Model
===================

    Implementation of a basic fragmentation model that allows for 
    correction of read depths to account for differences in mRNA
    length and edge effects due to sequencing being performed on 
    fragements of a certain average size and standard deviation.
    This is based on the approach taken by Ben Gordon, but rather
    than using a sampling approach attempts to calculate the 
    result analytically.

    Assumptions:

    1. Fragmentation is a uniformly random process along the mRNA.
          This is not necessarily valid for any fragmentation 
          method, but it holds well enough to account for the 
          majority of the variability. Additional baises can be
          added at a later date.

    2. Fragments are drawn from a Normal distribution defined by
          a mean length (bp) and standard deviation (bp). This is
          essential as it enables an exact solution to be found
          for the predicted factor at a given base in the sequence.
"""
#    Fragmentation Model
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import math
import csv
import random
import numpy as np
from scipy.stats import norm

def frag_factor_profile (mrna_len, frag_mean, frag_sd, bp_cutoff=None, no_end=False):
	"""Generate the expected fragmentation correction factor profile
	for a given length mRNA and fragments of given mean length and
	standard deviation.

	Parameters
	----------
	mrna_len : int
	    Length of the mRNA (bp).

	frag_mean : float
		Mean fragment length used during sequencing

	frag_sd : float
		Standard deviation in the fragment lengths.

	bp_cutoff : int (default=None)
		Limit of size of potential fragment. None assumes up to length of TU.

	no_end : bool (default=False)
		Flag whether the mRNA has an end.

	Returns
	-------
	profile: array
	    Profile of fragmentation factor along the mRNA.
	"""
	d = np.zeros(mrna_len)
	for bp in range(mrna_len):
		# Correct the potential mRNA lengths if no terminator
		corrected_mrna_len = mrna_len
		if no_end == True:
			corrected_mrna_len = mrna_len + frag_mean #+ bp_cutoff
		d[bp] = frag_at_base(bp, corrected_mrna_len, frag_mean, frag_sd, 
			                 bp_cutoff=bp_cutoff, no_end=no_end)
	return d

def frag_at_base (bp, mrna_len, frag_mean, frag_sd, bp_cutoff=None, no_end=False):
	"""Calculate the fragmentation factor at a given base in an mRNA.

	Parameters
	----------
	bp : int
		Base pair being considered in mRNA.

	mrna_len : int
	    Length of the mRNA (bp).

	frag_mean : float
		Mean fragment length used during sequencing

	frag_sd : float
		Standard deviation in the fragment lengths.

	bp_cutoff : int (default=None)
		Limit of size of potential fragment. None assumes up to length of TU.

	no_end : bool (default=False)
		Flag whether the mRNA has an end.

	Returns
	-------
	factor: float
	    Fragmentation factor at that base.
	"""
	result = 0.0
	min_bp = 1
	if bp_cutoff != None and bp_cutoff < frag_mean:
		min_bp = frag_mean - bp_cutoff
	max_bp = mrna_len+1
	if bp_cutoff != None and frag_mean+bp_cutoff <= mrna_len+1:
		max_bp = frag_mean + bp_cutoff
	for i in range(min_bp, max_bp):
		frag_cs = 0.0
		if no_end == False:
			frag_cs = frag_combinations(bp,i,mrna_len)
		else:
			frag_cs = frag_combinations(bp,i,mrna_len+max_bp)
		result += ( (1.0/(frag_sd*math.sqrt(2.0*math.pi))) * math.exp(
			       -(math.pow(i-2.0*frag_mean,2.0)
			       / (2.0*math.pow(frag_sd,2.0)))) 
		           * frag_cs )
	return result

def frag_c (i, x):
	"""Helper function to handle edge cases.

	Parameters
	----------
	i : int
	    Index of base in mRNA.

	x : int
		Fragment length being considered

	Returns
	-------
	combinations: int
	    Number of ways base can be positioned on one side.
	"""
	if i <= x:
		return 0.0
	else:
		return i-x

def frag_combinations (i, x, x_max):
	"""Number of ways a fragment of a given length can be uniquely positioned
	across a particular base in the mRNA.

	Parameters
	----------
	i : int
	    Index of base in mRNA.

	x : int
		Fragment length being considered

	x_max : int
		Maximum fragment length (length of mRNA).

	Returns
	-------
	combinations: int
	    Number of unique ways a fragment of that length can take across
	    that base.
	"""
	return i - frag_c(i,x) - frag_c(i,x_max-x)

def variant_frag_factor_profile (gcl, variant, frag_mean=280, frag_sd=70, bp_cutoff=None):
	"""Generate the expected fragmentation correction factor profile for all TUs in variant.

	Parameters
	----------
	gcl : GeneClusterLibrary
		Gene cluster library to use.

	variant : string
		Vairant name

	mrna_len : int
	    Length of the mRNA (bp).

	frag_mean : float
		Mean fragment length used during sequencing

	frag_sd : float
		Standard deviation in the fragment lengths.

	bp_cutoff : int (default=None)
		Limit of size of potential fragment. None assumes up to length of TU.

	Returns
	-------
	fwd_profiles: list
	    List containing elements of form [p_idx, t_idx, p_bp, t_bp, profile] for each
	    TU on the forward (+) strand.

	rev_profiles: list
	    List containing elements of form [p_idx, t_idx, p_bp, t_bp, profile] for each
	    TU on the reverse (-) strand.
	"""
	# Get the variant data
	var_data = gcl.variant_data(variant)
	# Extract all transcriptional units from the GeneClusterLibrary
	tus = gcl.transcriptional_units(non_terminated=True)
	# For each TU generate a profile
	fwd_profiles = []
	rev_profiles = []
	for tu in tus[variant]:
		p_bp = gcl.variant_part_idx_end_bp(variant, tu[0])
		tu_has_no_end = False
		t_bp = 0
		if tu[1] == None:
			tu_has_no_end = True
			if gcl.variant_part_idx_dir(variant, tu[0]) == 'F':
				tu[1] = len(var_data['part_list'])-1
				t_bp = gcl.variant_part_idx_end_bp(variant, tu[1])
			else:
				tu[1] = 0
				t_bp = gcl.variant_part_idx_start_bp(variant, tu[1])
		else:
			t_bp = gcl.variant_part_idx_start_bp(variant, tu[1])
		# Check direction and generate profile
		if p_bp < t_bp:
			profile = frag_factor_profile(t_bp-p_bp, frag_mean, frag_sd, 
				                          bp_cutoff=bp_cutoff, no_end=tu_has_no_end)
			# Forward strand
			fwd_profiles.append([tu[0], tu[1], p_bp, t_bp, profile])
		else:
			profile = frag_factor_profile(p_bp-t_bp, frag_mean, frag_sd,
				                          bp_cutoff=bp_cutoff, no_end=tu_has_no_end)
			# Forward strand
			rev_profiles.append([tu[0], tu[1], p_bp, t_bp, profile])
	return fwd_profiles, rev_profiles

def frag_profile_library (min_len, max_len, frag_mean=280, frag_sd=70, bp_cutoff=None):
	"""Generate library of fragmentation profiles over a specified range.

	Parameters
	----------
	min_len : int
	    Minimum length of an mRNA (bp).

	max_len : int
	    Maximum length of an mRNA (bp).

	frag_mean : float (default=280)
		Mean fragment length used during sequencing

	frag_sd : float (default=70)
		Standard deviation in the fragment lengths.

	bp_cutoff : int (default=None)
		Limit of size of potential fragment. None assumes up to length of TU.\

	Returns
	-------
	profiles_fixed: dict(list(float))
	    A dictionary of lists (the profiles) with the key corresponding to the length
	    of the mRNA. This dictionary considers mRNAs where the end it fixed.

	profiles_open: dict(list(float))
	    A dictionary of lists (the profiles) with the key corresponding to the length
	    of the mRNA. This dictionary considers mRNAs where the end it not fixed i.e.,
	    translation can continue outside of the constructs design. 
	"""
	profiles_fixed = {}
	profiles_open = {}
	for mrna_len in range(min_len, max_len+1):
		profiles_fixed[mrna_len] = frag_factor_profile(mrna_len, frag_mean, frag_sd, bp_cutoff=bp_cutoff, no_end=False)
		profiles_open[mrna_len] = frag_factor_profile(mrna_len, frag_mean, frag_sd, bp_cutoff=bp_cutoff, no_end=True)
	return profiles_fixed, profiles_open

def save_frag_profile_library (out_filename, min_len, max_len, frag_mean=280, frag_sd=70, bp_cutoff=None):
	"""Save a library of fragmentation profiles to a text file.

	Parameters
	----------
	out_filename : string
		Output filename

	min_len : int
	    Minimum length of an mRNA (bp).

	max_len : int
	    Maximum length of an mRNA (bp).

	frag_mean : float (default=280)
		Mean fragment length used during sequencing

	frag_sd : float (default=70)
		Standard deviation in the fragment lengths.

	bp_cutoff : int (default=None)
		Limit of size of potential fragment. None assumes up to length of TU.
	"""
	# Generate the profile library
	profiles_fixed, profiles_open = frag_profile_library(min_len, max_len, frag_mean=frag_mean, frag_sd=frag_sd, bp_cutoff=bp_cutoff)
	# Save the profiles to file
	f_out = open(out_filename, 'w')
	f_out.write('type,length,profile\n')
	for mrna_len in sorted(profiles_fixed.keys()):
		f_out.write('F,' + str(mrna_len) + ',' + ';'.join([str(x) for x in profiles_fixed[mrna_len]]) + '\n')
		f_out.write('O,' + str(mrna_len) + ',' + ';'.join([str(x) for x in profiles_open[mrna_len]]) + '\n')
	f_out.close()

def load_frag_profile_library (in_filename):
	"""Load a library of fragmentation profiles from a text file.

	Parameters
	----------
	in_filename : string
		Input filename

	Returns
	-------
	profiles_fixed: dict(list(float))
	    A dictionary of lists (the profiles) with the key corresponding to the length
	    of the mRNA. This dictionary considers mRNAs where the end it fixed.

	profiles_open: dict(list(float))
	    A dictionary of lists (the profiles) with the key corresponding to the length
	    of the mRNA. This dictionary considers mRNAs where the end it not fixed i.e.,
	    translation can continue outside of the constructs design. 
	"""
	profiles_fixed = {}
	profiles_open = {}
	# Load the file data
	file_reader = csv.reader(open(in_filename, 'rU'), delimiter=',')
	# Ignore the header
	header = next(file_reader)
	# Load the rest of the file
	for row in file_reader:
		if row[0] == 'F':
			profiles_fixed[int(row[1])] = [float(x) for x in row[2].split(';')]
		else:
			profiles_open[int(row[1])] = [float(x) for x in row[2].split(';')]
	return profiles_fixed, profiles_open

def frag_factor_profile_random_prob_model (mrna_len, frag_mean, frag_sd, mrna_count=10000, bp_frag_prob=0.01, verbose=False):
	# The list of fragments we generate (tuples of start and end bp)
	frags = []
	# Randomly fragment the mRNAs with uniform probability and probabistically select those of the right size
	for mrna in range(mrna_count):
		if verbose==True:
			print 'Processing mRNA:', mrna
		start_bp = 0
		end_bp = 0
		for bp in range(mrna_len):
			if random.random() <= bp_frag_prob:
				# Fragment is made
				end_bp = bp+1
				frag_len = end_bp-start_bp
				frag_prob = norm.pdf(frag_len, loc=frag_mean, scale=frag_sd)
				if random.random() <= frag_prob:
					# The fragment is selected
					frags.append([start_bp, end_bp])
				start_bp=end_bp
		if end_bp < mrna_len:
			# Add the last fragment
			frag_len = mrna_len-start_bp
			frag_prob = norm.pdf(frag_len, loc=frag_mean, scale=frag_sd)
			if random.random() <= frag_prob:
				# The fragment is selected
				frags.append([start_bp, mrna_len])
	if verbose==True:
		print 'Found', len(frags), 'fragments.'
	# Count up fragments that fall at each point along mRNA
	frag_profile = np.zeros(mrna_len)
	for frag in frags:
		c = np.zeros(mrna_len)
		c[frag[0]:frag[1]] = np.ones(frag[1]-frag[0])
		frag_profile = frag_profile + c
	return frag_profile

def frag_factor_profile_random_break_model (mrna_len, frag_mean, frag_sd, mrna_count=10000, max_frags=10, verbose=False):
	# The list of fragments we generate (tuples of start and end bp)
	frags = []
	# Randomly fragment the mRNAs with uniform probability and probabistically select those of the right size
	for mrna in range(mrna_count):
		if verbose==True:
			print 'Processing mRNA:', mrna
		# Break the mRNA into pieces
		breaks = []
		num_of_breaks = random.randint(1,max_frags+1)
		for b in range(num_of_breaks):
			b_bp = random.randint(0,mrna_len)
			while b_bp in breaks:
				b_bp = random.randint(0,mrna_len)
			breaks.append(b_bp)

		# Probabilistically choose fragments based on distirbution
		start_bp = 0
		for b in breaks:
			end_bp = b+1
			frag_len = end_bp-start_bp
			frag_prob = norm.pdf(frag_len, loc=frag_mean, scale=frag_sd)
			if random.random() <= frag_prob:
				# The fragment is selected
				frags.append([start_bp, end_bp])
			start_bp=end_bp
			if end_bp < mrna_len:
				# Add the last fragment
				frag_len = mrna_len-start_bp
				frag_prob = norm.pdf(frag_len, loc=frag_mean, scale=frag_sd)
				if random.random() <= frag_prob:
					# The fragment is selected
					frags.append([start_bp, mrna_len])
	if verbose==True:
		print 'Found', len(frags), 'fragments.'
	# Count up fragments that fall at each point along mRNA
	frag_profile = np.zeros(mrna_len)
	for frag in frags:
		c = np.zeros(mrna_len)
		c[frag[0]:frag[1]] = np.ones(frag[1]-frag[0])
		frag_profile = frag_profile + c
	return frag_profile


