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

def scar_compatible (scar_set, new_scar, max_homology=3, allowed_set=[]):
	"""Check if scar compaible with current set at given level of homology.
	"""
	# Check if not allowed... then test properties
	if new_scar not in allowed_set:
		# Check not all a single base
		scar_len = len(new_scar)
		if new_scar == 'A'*scar_len or new_scar == 'T'*scar_len or new_scar == 'G'*scar_len or new_scar == 'C'*scar_len:
			return False
		# Check G/C content
		g_count = new_scar.count('G')
		c_count = new_scar.count('C')
		gc_per = float(g_count+c_count)/len(new_scar)
		if gc_per < 0.25 or gc_per > 0.75:
			return False
	if scar_set != []:
		# Must not also be in the reverse complement set so generate for checking
		full_check_set = list(scar_set)
		for scar in scar_set:
			full_check_set.append(reverse_complement(scar))
		for scar in full_check_set:
			homology_count = 0
			for scar_bp in range(len(scar)):
				if scar[scar_bp] == new_scar[scar_bp]:
					homology_count += 1
			if homology_count > max_homology:
				return False
	return True

# Complement function: http://stackoverflow.com/questions/19570800/reverse-complement-dna
def complement(s): 
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	letters = list(s) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters)

# Reverse complement function: http://stackoverflow.com/questions/19570800/reverse-complement-dna
def reverse_complement(s):
	return complement(s[::-1])

def random_scar (scar_len, exclude_set=[], max_homology=3, allowed_set=[]):
	"""Generate a random scar of a given length that is not in the excluded set.
	"""
	# Used to ensure we don't get involved in infinite loop
	i = 0
	# Must not also be in the reverse complement set so generate for checking
	rc_exclude_set = []
	for scar in exclude_set:
		rc_exclude_set.append(reverse_complement(scar))
	# Randomly generate candidates until valid scar found (or maximum iterations reached)
	while True:
		# From: http://stackoverflow.com/questions/21205836/generating-random-sequences-of-dna
		new_scar = ''.join(random.choice('ATGC') for _ in xrange(scar_len))
		if scar_compatible(exclude_set, new_scar, max_homology=max_homology, allowed_set=allowed_set) == False:
			break
		i += 1
		if i > MAX_RANDOM_ITERATIONS:
			new_scar = None
			break
	return new_scar

def enumerate_scars (scar_len, depth, cur_scar, full_scar_set, new_scars, max_homology=3, found=[0], num_to_find=None, allowed_set=[]):
	"""Recursive method to enumerate all possible scars (be careful memory usage high for large lengths)
	"""
	if depth == scar_len:
		if num_to_find == None or found[0] < num_to_find:
			if scar_compatible(full_scar_set, cur_scar, max_homology=max_homology, allowed_set=allowed_set) == True:
					full_scar_set.append(cur_scar)
					new_scars.append(cur_scar)
					found[0] = found[0] + 1
	else:
		bases = list('ATGC')
		for base in bases:
			enumerate_scars(scar_len, depth + 1, cur_scar + base, full_scar_set, new_scars, max_homology=max_homology, found=found, num_to_find=num_to_find, allowed_set=allowed_set)

def find_scars (scar_len, seed_set=[], max_homology=3, num_to_find=None, random_search=True, allowed_set=[]):
	"""Find a required number of scars orthogonal for a seed set with a 
	maximim number of bp homology.
	"""
	# 1. Check pairwise all scars in seed set and remove if issues (track removed)
	seed_set = [x.upper() for x in seed_set] # We only work in capitals
	seed_set = list(set(seed_set)) # Remove duplicates
	# We only work with upper case sequences
	allowed_set = [x.upper() for x in allowed_set]
	removed_from_seed = []
	full_scar_set = []
	new_scars = []
	for cur_scar in seed_set:
		if scar_compatible(full_scar_set, cur_scar, max_homology=max_homology, allowed_set=allowed_set) == True:
			full_scar_set.append(cur_scar)
		else:
			removed_from_seed.append(cur_scar)
	# 2. Search for orthogonal scars (use random search or enumerate all)
	if num_to_find == None and random_search == True:
		# Only perform random search for less than all scars
		random_search = False
	if random_search == True:
		# Perform a random search (gives greater variability in output set)
		for i in range(num_to_find):
			new_scar = random_scar(scar_len, full_scar_set, max_homology=max_homology, allowed_set=allowed_set)
			if new_scar == None:
				# Could not find scar exist
				return None, None, None
			else:
				full_scar_set.append(new_scar)
				new_scars.append(new_scar)
	else:
		# Enumerate all possibilities until all required number of scars are found (or all)
		if num_to_find == None:
			enumerate_scars(scar_len, 0, '', full_scar_set, new_scars, max_homology=max_homology, allowed_set=allowed_set)
		else:
			found = [0]
			enumerate_scars(scar_len, 0, '', full_scar_set, new_scars, max_homology=max_homology, found=found, num_to_find=num_to_find, allowed_set=allowed_set)
			if found[0] < num_to_find:
				return None, None, None
	return full_scar_set, new_scars, removed_from_seed
