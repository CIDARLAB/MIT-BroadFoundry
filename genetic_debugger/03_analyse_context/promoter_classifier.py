#!/usr/bin/env python
"""Promoter classifier
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import gene_cluster_library as gcl
import numpy as np
import csv

##################################################################################
# PROMOTER CONTEXT CLASSIFICATION
##################################################################################

class PromoterClassifier ():

	# Single promoter contexts with classify
	PROMOTER_TYPES = ['P-CDS', 'T-P-CDS', 'CDS-P-CDS',
	                  'P-P', 'P-P-P', 'CDS-P-P', 'T-P-P', 'P-P-CDS'
	                  'rP-P-CDS']

	# Major context types
	SINGLE_PROMOTER     = ['P-CDS', 'T-P-CDS', 'CDS-P-CDS']
	MULTI_PROMOTER      = ['P-P', 'P-P-P', 'CDS-P-P', 'T-P-P', 'P-P-CDS']
	DIVERGENT_PROMOTER  = ['rP-P-CDS']

	def classify (self, gcl):
		"""Run the clasification.
		"""
		# Somewhere to hold the classified promoter data
		p_data = []
		p_bad_data = []
		# Find all promoters in the library
		p_insts = gcl.find_part_type_instances('Promoter')
		for v_key in p_insts:
			for p_inst in p_insts[v_key]:
				# Collect some data about the promoter
				p_name = gcl.variant_part_idx_name(v_key, p_inst)
				p_dir = gcl.variant_part_idx_dir(v_key, p_inst)
				p_type = 'P'
				p_details = {}
				us_CDS_idx = None
				ds_CDS_idx = None
				us_Promoter_idx = None
				ds_Promoter_idx = None
				us_Terminator_idx = None
				us_CDS_name = None
				ds_CDS_name = None
				us_Promoter_name = None
				ds_Promoter_name = None
				us_Terminator_name = None
				p_include = True
				# Clasify by moving upstream and see what we find
				i = 0
				if p_dir == 'F':
					i = p_inst-1
				else:
					i = p_inst+1
				while i >= 0 and i < len(gcl.variant_part_list(v_key)):
					if gcl.variant_part_idx_type(v_key, i) == 'Promoter':
						if gcl.variant_part_idx_dir(v_key, i) == p_dir:
							p_type = 'P-' + p_type
							us_Promoter_idx = i
							us_Promoter_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							p_type = 'rP-' + p_type
							us_Promoter_idx = i
							us_Promoter_name = gcl.variant_part_idx_name(v_key, i)
							break
					if gcl.variant_part_idx_type(v_key, i) == 'Terminator':
						if gcl.variant_part_idx_dir(v_key, i) == p_dir:
							p_type = 'T-' + p_type
							us_Terminator_idx = i
							us_Terminator_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							p_include = False
							break
					if gcl.variant_part_idx_type(v_key, i) == 'CDS':
						if gcl.variant_part_idx_dir(v_key, i) == p_dir:
							p_type = 'CDS-' + p_type
							us_CDS_idx = i
							us_CDS_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							p_include = False
							break
					# Move to next (previous) element
					if p_dir == 'F':
						i -= 1
					else:
						i += 1

				# Classify downstream by moving and see what we get
				if p_dir == 'F':
					i = p_inst+1
				else:
					i = p_inst-1
				while i >= 0 and i < len(gcl.variant_part_list(v_key)):
					if gcl.variant_part_idx_type(v_key, i) == 'Promoter':
						if gcl.variant_part_idx_dir(v_key, i) == p_dir:
							p_type = p_type + '-P'
							ds_Promoter_idx = i
							ds_Promoter_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							p_include = False
							break
					if gcl.variant_part_idx_type(v_key, i) == 'Terminator':
						p_include = False
						break
					if gcl.variant_part_idx_type(v_key, i) == 'CDS':
						if gcl.variant_part_idx_dir(v_key, i) == p_dir:
							p_type = p_type + '-CDS'
							ds_CDS_idx = i
							ds_CDS_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							p_include = False
							break
					# Move to next (previous) element
					if p_dir == 'F':
						i += 1
					else:
						i -= 1
				# Check that valid type was found
				p_details['us_CDS_idx'] = us_CDS_idx
				p_details['ds_CDS_idx'] = ds_CDS_idx
				p_details['us_Promoter_idx'] = us_Promoter_idx
				p_details['ds_Promoter_idx'] = ds_Promoter_idx
				p_details['us_Terminator_idx'] = us_Terminator_idx
				p_details['us_CDS_name'] = us_CDS_name
				p_details['ds_CDS_name'] = ds_CDS_name
				p_details['us_Promoter_name'] = us_Promoter_name
				p_details['ds_Promoter_name'] = ds_Promoter_name
				p_details['us_Terminator_name'] = us_Terminator_name
				if p_include == True:
					p_data.append([v_key, p_inst, p_name, p_type, p_details])
				else:
					p_bad_data.append([v_key, p_inst, p_name, p_type, p_details])
		return p_data, p_bad_data

##################################################################################
# PROMOTER STRENGTH CALCULATIONS
##################################################################################

class PromoterStrengthFromTrace ():

	def __init__ (self, fwd_skip=0, fwd_len=30, rev_skip=0, rev_len=30):
		self.fwd_skip = fwd_skip
		self.fwd_len = fwd_len
		self.rev_skip = rev_skip
		self.rev_len = rev_len

	def append_strength (self, gcl, traces, p_data):
		# For each promoter calcuate the strength and append to p_data element
		for p_inst in p_data:
			v_key = p_inst[0]
			p_idx = p_inst[1]
			p_details = p_inst[4]
			trace = traces[v_key]
			p_strength = self.promoter_strength(gcl, v_key, p_idx, 
				                                p_details['ds_CDS_idx'], 
				                                p_details['us_CDS_idx'],
				                                traces[v_key])
			p_inst.append(p_strength)

	def promoter_strength (self, gcl, variant, p_idx, ds_CDS_idx, us_CDS_idx, trace):
		"""Calculate an absolute promoter strength based on transcriptomic data.
		"""
		var_data = gcl.variant_data(variant)
		pre_tx = []
		post_tx = []
		# Find the direction of the TU
		step = 1
		if gcl.variant_part_idx_dir(variant, p_idx) != 'F':
			step = -1
		p_start_idx = gcl.variant_part_idx_seq_idx(variant, p_idx)
		# Strand alters where data comes from, extract relevant information
		if step == 1:
			# +ve strand
			p_end_idx = p_start_idx + gcl.variant_part_idx_seq_len(variant, p_idx)
			pre_tx = trace[0][p_start_idx-self.rev_skip-self.rev_len:p_start_idx-self.rev_skip]
			post_tx = trace[0][p_end_idx+self.fwd_skip:p_end_idx+self.fwd_skip+self.fwd_len]
		else:
			# -ve strand
			p_end_idx = p_start_idx - gcl.variant_part_idx_seq_len(variant, p_idx)
			pre_tx = trace[1][p_start_idx+self.rev_skip:p_start_idx+self.rev_skip+self.rev_len]
			post_tx = trace[1][p_end_idx-self.fwd_skip-self.fwd_len:p_end_idx-self.fwd_skip]
		# Promoter strength is the abs change in read depth (median used to smooth noise)
		p_strength = np.median(post_tx) - np.median(pre_tx)
		return [np.median(pre_tx), np.median(post_tx), p_strength]

class PromoterStrengthFromRSEM ():

	def __init__ (self, rsem_filename):
		rsem_data = {}
		file_reader = csv.reader(open(rsem_filename, 'rU'), delimiter=',')
		# Process header for column names (used in dict returned)
		header = next(file_reader)
		# Process each line
		for row in file_reader:
			if len(row) == 5:
				p_variant = row[0]
				p_idx = int(row[1])
				p_name = row[2]
				p_tx_in = float(row[3])
				p_tx_out = float(row[4])
				if p_name not in rsem_data.keys():
					rsem_data[p_name] = []
				rsem_data[p_name].append([p_variant, p_idx, [p_tx_in,p_tx_out], p_tx_out])
		self.rsem_data = rsem_data

	def append_strength (self, p_data):
		# For each terminator calcuate the strength and append to t_data element
		for p_inst in p_data:
			v_key = p_inst[0]
			p_idx = p_inst[1]
			p_name = p_inst[2]
			p_strength = self.promoter_strength(v_key, p_idx, p_name)
			p_inst.append(p_strength)

	def promoter_strength (self, variant, p_idx, p_name):
		"""Calculate an absolute promoter strength based on RSEM data.
		"""
		found_data = [-1.0, -1.0, -1.0]
		for p_inst in self.rsem_data[p_name]:
			if p_inst[0] == variant and p_inst[1] == p_idx:
				found_data = [p_inst[2][0], p_inst[2][1], p_inst[3]]
				break
		return found_data
