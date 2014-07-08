#!/usr/bin/env python
"""
Scar Choice
===========

    Script to help with the selection of orthogonal scars used for 
    Golden-Gate like assembly of constructs.
"""
#    Scar Choice
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import random

MAX_RANDOM_ITERATIONS = 10000

def scar_compatible (scar_set, new_scar, max_homology=3):
	"""Check if scar compaible with current set at given level of homology.
	"""
	if scar_set != []:
		for scar in scar_set:
			homology_count = 0
			for scar_bp in range(len(scar)):
				if scar[scar_bp] == new_scar[scar_bp]:
					homology_count += 1
			if homology_count > max_homology:
				return False
	return True

def random_scar (scar_len, exclude_set=[]):
	"""Generate a random scar of a given length that is not in the excluded set.
	"""
	# From: http://stackoverflow.com/questions/21205836/generating-random-sequences-of-dna
	i = 0
	while true:
		new_scar = ''.join(random.choice('CGTA') for _ in xrange(scar_len))
		# Check it is not one in the excluded list
		if new_scar != 'A'*scar_len and new_scar != 'T'*scar_len and
		   new_scar != 'G'*scar_len and new_scar != 'C'*scar_len and 
		   new_scar not in exclude_set:
			break
		i += 1
		if i > MAX_RANDOM_ITERATIONS:
			new_scar = None
			break
	return new_scar

def find_scars (seed_set=[], scar_len, max_homology=3, num_to_find=None, random_search=True):
	"""Find a required number of scars orthogonal for a seed set with a 
	maximim number of bp homology.
	"""
	# Check that number to find is valid
	# TODO
	# 1. Check if seed needs to be generated
	if seed_set == []:
		new_scar = random_scar(scar_len)
		if new_scar == None:
			return None, None, None
		else:
			seed_set = [new_scar]
	# 2. Check pairwise all scars in seed set and remove if issues (track removed)
	seed_set = list(set(seed_set)) # Remove duplicates
	removed_from_seed = []
	full_scar_set = []
	for cur_scar in range(len(seed_set)):
		if scar_compatible(full_scar_set, cur_scar, max_homology=max_homology):
			full_set.append(cur_scar)
		else:
			removed_from_seed(cur_scar)
	# 3. Search for orthogonal scars (use random search or enumerate all)
	if random_search:
		# 

	else:
		# Enumerate all possibilities until number to found matched

	return full_scar_set, new_scars, removed_from_seed








