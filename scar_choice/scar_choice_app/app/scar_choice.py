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
	scar_len = len(new_scar)
	if new_scar == 'A'*scar_len or new_scar == 'T'*scar_len or new_scar == 'G'*scar_len or new_scar == 'C'*scar_len:
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

# Compliment function: http://stackoverflow.com/questions/19570800/reverse-complement-dna
def complement(s): 
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	letters = list(s) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters)

# Reverse compliment function: http://stackoverflow.com/questions/19570800/reverse-complement-dna
def reverse_complement(s):
	return complement(s[::-1])

def random_scar (scar_len, exclude_set=[]):
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
		# Check it is not one in the excluded list
		if new_scar != 'A'*scar_len and new_scar != 'T'*scar_len and new_scar != 'G'*scar_len and new_scar != 'C'*scar_len and new_scar not in exclude_set and new_scar not in rc_exclude_set:
			break
		i += 1
		if i > MAX_RANDOM_ITERATIONS:
			new_scar = None
			break
	return new_scar

def enumerate_scars (scar_len, depth, cur_scar, full_scar_set, new_scars, max_homology=3, found=[0], num_to_find=None):
	"""Recursive method to enumerate all possible scars (be careful memory usage high for large lengths)
	"""
	if depth == scar_len:
		if num_to_find == None or found[0] < num_to_find:
			if scar_compatible(full_scar_set, cur_scar, max_homology=max_homology):
					full_scar_set.append(cur_scar)
					new_scars.append(cur_scar)
					found[0] = found[0] + 1
	else:
		bases = list('ATGC')
		for base in bases:
			enumerate_scars(scar_len, depth + 1, cur_scar + base, full_scar_set, new_scars, max_homology=max_homology, found=found, num_to_find=num_to_find)

def find_scars (scar_len, seed_set=[], max_homology=3, num_to_find=None, random_search=True):
	"""Find a required number of scars orthogonal for a seed set with a 
	maximim number of bp homology.
	"""
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
	new_scars = []
	for cur_scar in seed_set:
		if scar_compatible(full_scar_set, cur_scar, max_homology=max_homology):
			full_scar_set.append(cur_scar)
		else:
			removed_from_seed(cur_scar)
	# 3. Search for orthogonal scars (use random search or enumerate all)
	if num_to_find == None and random_search:
		# Only perform random search for less than all scars
		random_search = False
	if random_search:
		# Perform a random search (gives greater variability in output set)
		for i in range(num_to_find):
			new_scar = random_scar(scar_len, full_scar_set)
			if new_scar == None:
				# Could not find scar exist
				return None, None, None
			else:
				full_scar_set.append(new_scar)
				new_scars.append(new_scar)
	else:
		# Enumerate all possibilities until all required number of scars are found (or all)
		cur_scar = seed_set[0]
		new_scar = cur_scar
		if num_to_find == None:
			enumerate_scars(scar_len, 0, '', full_scar_set, new_scars, max_homology=max_homology)
		else:
			found = [0]
			enumerate_scars(scar_len, 0, '', full_scar_set, new_scars, max_homology=max_homology, found=found, num_to_find=num_to_find)
			if found[0] < num_to_find:
				return None, None, None
	return full_scar_set, new_scars, removed_from_seed
