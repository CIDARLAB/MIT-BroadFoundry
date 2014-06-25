#!/usr/bin/env python
"""
Gene Cluster Library Analysis Tools
===================================

    Gene Cluster Analysis Tools are a set of methods to robustly perform 
    general queries of gene cluster libraries.
"""
#    Gene Cluster Tools for Python
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import sys
if sys.version_info[:2] < (2, 6):
    m = "Python version 2.6 or later is required for Gene Cluster Tools (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

import csv
import gene_cluster_library as gcl

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def abs_promoter_strength (gcl, variant, p_idx, t_idx, tx_profile):
	"""Calculate an absolute promoter strength based on transcriptomic data.
	"""
	return None

def abs_promoter_strengths (gcl, tu_insts, tx_profiles):
	"""Calculate an absolute promoter strengths of a set of transcriptional units.
	"""
	return None

def filter_tu_insts (gcl, tu_insts, filter_insts, idx=0):
	results = {}
	for v_key in tu_insts.keys():
		if v_key in filter_insts.keys():
			new_tu_data = []
			tu_data = tu_insts[v_key]
			filter_data = filter_insts[v_key]
			for tu_val in tu_data:
				filter_flag = False
				for f_val in filter_data:
					if tu_val[idx] == f_val[0] or tu_val[idx] == f_val[1]:
						# This item should be filtered
						filter_flag = True
						break
				if filter_flag == False:
					new_tu_data.append(list(tu_val))
			results[v_key] = new_tu_data
		else:
			results[v_key] = list(tu_insts[v_key])
	return results

def filter_tu_insts_on_promoters (gcl, tu_insts, filter_insts):
	return filter_tu_insts(gcl, tu_insts, filter_insts, idx=0)

def filter_tu_insts_on_terminators (gcl, tu_insts, filter_insts):
	return filter_tu_insts(gcl, tu_insts, filter_insts, idx=1)

def tu_meta_data (gcl, variant, tu):
	"""Find meta data for the transcriptional unit.

	If any of the following components are not found then None is returned.

	Parameters
	----------
	gcl : GeneClusterLibrary
		Gene cluster library to perform query using.

	variant : string
		Gene cluster variant name.

	tu : list([start_part_idx, end_part_idx])
		Transcriptional unit start and end part indexes.

	Returns
	-------
	metadata : list(string) ([promoter, rbs, cds, terminator])
		Metadata for the transcriptional unit (part names)
	"""
	p = None
	r = None
	c = None
	t = None
	# Make sure we cycle in the right direction
	step_dir = 1
	if tu[1] < tu[0]:
	    step_dir = -1
	# Cycle through all parts to extract meta data
	for cur_part_idx in range(tu[0], tu[1]+step_dir, step_dir):
		cur_name = gcl.variants[variant]['part_list'][cur_part_idx]['part_name']
		cur_type = gcl.parts[cur_name]['type']
		# Check to see if found relevant part
		if p == None and cur_type == 'Promoter':
			p = cur_name
		if r == None and cur_type == 'RBS':
			r = cur_name
		if c == None and cur_type == 'CDS':
			c = cur_name
		if t == None and cur_type == 'Terminator':
			t = cur_name
	return [p, r, c, t]

def tu_insts_meta_data (gcl, tu_insts):
	"""Find meta data for a set of transcriptional unit instances.

	If any of the following components are not found then None is returned.

	Parameters
	----------
	gcl : GeneClusterLibrary
		Gene cluster library to perform query using.

	variant : string
		Gene cluster variant name.

	tu : list([start_part_idx, end_part_idx])
		Transcriptional unit start and end part indexes.

	Returns
	-------
	metadata : dict(list(string)) ([promoter, rbs, cds, terminator])
		Metadata for the transcriptional units (part names) with dictionary
		keyed by variant name.
	"""
	results = {}
	for v_key in tu_insts.keys():
		results[v_key] = []
		for tu in tu_insts[v_key]:
			results[v_key].append(tu_meta_data(gcl, v_key, tu))
	return results
