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
		for scar in scar_set:
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
	# From: http://stackoverflow.com/questions/21205836/generating-random-sequences-of-dna
	i = 0
	# Must not also be in the reverse complement set so generate for cheking
	rc_exclude_set = []
	for scar in exclude_set:
		rc_exclude_set.append(reverse_complement(scar))
	# Randomly generate candidates until valid scar found (or maximum iterations reached)
	while True:
		new_scar = ''.join(random.choice('ATGC') for _ in xrange(scar_len))
		# Check it is not one in the excluded list
		if new_scar != 'A'*scar_len and new_scar != 'T'*scar_len and new_scar != 'G'*scar_len and new_scar != 'C'*scar_len and new_scar not in exclude_set and new_scar not in rc_exclude_set:
			break
		i += 1
		if i > MAX_RANDOM_ITERATIONS:
			new_scar = None
			break
	return new_scar

def find_scars (scar_len, seed_set=[], max_homology=3, num_to_find=None, random_search=True):
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
	new_scars = []
	for cur_scar in seed_set:
		if scar_compatible(full_scar_set, cur_scar, max_homology=max_homology):
			full_scar_set.append(cur_scar)
		else:
			removed_from_seed(cur_scar)
	# 3. Search for orthogonal scars (use random search or enumerate all)
	if num_to_find == None and random_search:
		# Only perform random search for less that all scars
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
		# Enumerate all possibilities until all required number of scars are found
		cur_scar = seed_set[0]
		new_scar = cur_scar
		scars_found = 0
		for cur_idx in range(len(cur_scar)):
			for base in 'ATGC':
				temp = list(new_scar)
				temp[cur_idx] = base
				new_scar = ''.join(temp)
				if base != cur_scar[cur_idx] and scar_compatible(full_scar_set, new_scar, max_homology=max_homology):
					full_scar_set.append(new_scar)
					new_scars.append(new_scar)
					scars_found += 1
					if num_to_find != None and scars_found == num_to_find:
						return full_scar_set, new_scars, removed_from_seed
		if num_to_find != None and scars_found < num_to_find:
			return None, None, None
	return full_scar_set, new_scars, removed_from_seed
