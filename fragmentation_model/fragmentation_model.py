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
          essential as it enables and exact solution to be found
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
import numpy as np

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
			corrected_mrna_len = mrna_len + frag_mean + bp_cutoff
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
