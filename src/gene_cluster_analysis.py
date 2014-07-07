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
import numpy as np
import gene_cluster_library as gcl

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def abs_promoter_strength (gcl, variant, p_idx, t_idx, tx_profile, fwd_skip=200, 
	                       fwd_len=100, rev_skip=0, rev_len=30):
	"""Calculate an absolute promoter strength based on transcriptomic data.
	"""
	var_data = gcl.variants[variant]
	pre_tx = []
	post_tx = []
	# Find the direction of the TU
	step = 1
	if t_idx < p_idx:
		step = -1
	p_start_idx = var_data['part_list'][p_idx]['seq_idx']
	# Strand alters where data comes from, extract relevant information
	if step == 1:
		# +ve strand
		p_end_idx = p_start_idx + var_data['part_list'][p_idx]['seq_len']
		pre_tx = tx_profile[0][p_start_idx-rev_skip-rev_len:p_start_idx-rev_skip]
		post_tx = tx_profile[0][p_end_idx+fwd_skip:p_end_idx+fwd_skip+fwd_len]
	else:
		# -ve strand
		p_end_idx = p_start_idx - var_data['part_list'][p_idx]['seq_len']
		pre_tx = tx_profile[1][p_start_idx+rev_skip:p_start_idx+rev_skip+rev_len]
		post_tx = tx_profile[1][p_end_idx-fwd_skip-fwd_len:p_end_idx-fwd_skip]
	# Promoter strength is the abs change in read depth (median used to smooth noise)
	p_strength = np.median(post_tx) - np.median(pre_tx)
	return p_strength

def abs_promoter_strengths (gcl, tu_insts, tx_profiles, fwd_skip=200, fwd_len=100, 
	                        rev_skip=0, rev_len=30):
	"""Calculate an absolute promoter strengths for a set of transcriptional units.
	"""
	results = {}
	for v_key in tu_insts.keys():
		results[v_key] = []
		for tu in tu_insts[v_key]:
			# For each TU, calculate strength and store
			results[v_key].append(abs_promoter_strength(gcl, v_key, tu[0], tu[1], 
				                  tx_profiles[v_key], fwd_skip=fwd_skip, fwd_len=fwd_len, 
				                  rev_skip=rev_skip, rev_len=rev_len))
	return results

def attrib_promoter_strength (gcl, variant, p_idx, attrib_name='rpkm'):
	"""Calculate promoter strength based on attribute data.
	"""
	# Find the CDS downstream of the promoter
	cds_idx = gcl.find_next_part_idx(variant, p_idx, part_type='CDS')
	# Use the attrib_name attribute as strength
	p_strength = gcl.variants[variant]['part_list'][cds_idx][attrib_name]
	return p_strength

def attrib_promoter_strengths (gcl, tu_insts, attrib_name='rpkm'):
	"""Calculate promoter strength based on attribute data for a set of 
	transcriptional units.
	"""
	results = {}
	for v_key in tu_insts.keys():
		results[v_key] = []
		for tu in tu_insts[v_key]:
			# For each TU, calculate strength and store
			results[v_key].append(attrib_promoter_strength(gcl, v_key, tu[0],
				                                           attrib_name=attrib_name))
	return results

def filter_tu_insts (tu_insts, filter_insts, idx=0):
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

def filter_tu_insts_on_promoters (tu_insts, filter_insts):
	return filter_tu_insts(tu_insts, filter_insts, idx=0)

def filter_tu_insts_on_terminators (tu_insts, filter_insts):
	return filter_tu_insts(tu_insts, filter_insts, idx=1)

def load_strand_data (filename):
	data = {}
	reader = csv.reader(open(filename, 'rb'), delimiter=',')
	# Use the header to get the names of the keys for the prediction data
	header = reader.next()
	for el in header[1:]:
		el_name = el[6:-1]
		# Only consider the first of the replicates
		if el_name.count('.') == 1:
			el_name = el_name[0]
		if el_name not in data.keys():
			data[el_name] = [[],[]]
	for row in reader:
		for el_idx in range(1, len(row)):
			el_name = header[el_idx][6:-1]
			# Only consider the first of the replicates
			if el_name.count('.') == 1:
				el_name = el_name[0]
			strand = header[el_idx][-1]
			if strand == '+':
				data[el_name][0].append(float(row[el_idx]))
			else:
				data[el_name][1].append(float(row[el_idx]))
	return data
