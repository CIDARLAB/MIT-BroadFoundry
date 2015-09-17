#!/usr/bin/env python
"""Terminator classifier
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import gene_cluster_library as gcl
import numpy as np
import csv

##################################################################################
# TERMINATOR CONTEXT CLASSIFICATION
##################################################################################

class TerminatorClassifier ():

	# Single terminator contexts with classify
	TERMINATOR_TYPES = ['T-P', 'CDS-T-P', 'CDS-T-T', 'T-T-P',
	                    'CDS-T-rT', 'CDS-T-CDS', 'CDS-T-rCDS', 'T-rCDS']

	def classify (self, gcl):
		"""Run the clasification.
		"""
		# Somewhere to hold the classified terminator data
		t_data = []
		t_bad_data = []
		# Find all terminators in the library
		t_insts = gcl.find_part_type_instances('Terminator')
		for v_key in t_insts:
			for t_inst in t_insts[v_key]:
				# Collect some data about the terminator
				t_name = gcl.variant_part_idx_name(v_key, t_inst)
				t_dir = gcl.variant_part_idx_dir(v_key, t_inst)
				t_type = 'T'
				t_details = {}
				us_CDS_idx = None
				ds_CDS_idx = None
				us_Promoter_idx = None
				ds_Promoter_idx = None
				us_Terminator_idx = None
				ds_Terminator_idx = None
				us_CDS_name = None
				ds_CDS_name = None
				us_Promoter_name = None
				ds_Promoter_name = None
				us_Terminator_name = None
				ds_Terminator_name = None
				t_include = True
				# Clasify by moving upstream and see what we find
				i = 0
				if t_dir == 'F':
					i = t_inst-1
				else:
					i = t_inst+1
				while i >= 0 and i < len(gcl.variant_part_list(v_key)):
					if gcl.variant_part_idx_type(v_key, i) == 'Promoter':
						us_Promoter_idx = i
						us_Promoter_name = gcl.variant_part_idx_name(v_key, i)
						t_include = False
						break
					if gcl.variant_part_idx_type(v_key, i) == 'Terminator':
						if gcl.variant_part_idx_dir(v_key, i) == t_dir:
							t_type = 'T-' + t_type
							us_Terminator_idx = i
							us_Terminator_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							us_Terminator_idx = i
							us_Terminator_name = gcl.variant_part_idx_name(v_key, i)
							t_include = False
							break
					if gcl.variant_part_idx_type(v_key, i) == 'CDS':
						if gcl.variant_part_idx_dir(v_key, i) == t_dir:
							t_type = 'CDS-' + t_type
							us_CDS_idx = i
							us_CDS_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							us_CDS_idx = i
							us_CDS_name = gcl.variant_part_idx_name(v_key, i)
							t_include = False
							break
					# Move to next (previous) element
					if t_dir == 'F':
						i -= 1
					else:
						i += 1

				# Classify downstream by moving and see what we get
				if t_dir == 'F':
					i = t_inst+1
				else:
					i = t_inst-1
				while i >= 0 and i < len(gcl.variant_part_list(v_key)):
					if gcl.variant_part_idx_type(v_key, i) == 'Promoter':
						if gcl.variant_part_idx_dir(v_key, i) == t_dir:
							t_type = t_type + '-P'
							ds_Promoter_idx = i
							ds_Promoter_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							t_type = t_type + '-rP'
							ds_Promoter_idx = i
							ds_Promoter_name = gcl.variant_part_idx_name(v_key, i)
							break
					if gcl.variant_part_idx_type(v_key, i) == 'Terminator':
						if gcl.variant_part_idx_dir(v_key, i) == t_dir:
							t_type = t_type + '-T'
							ds_Terminator_idx = i
							ds_Terminator_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							t_type = t_type + '-rT'
							ds_Terminator_idx = i
							ds_Terminator_name = gcl.variant_part_idx_name(v_key, i)
							break
					if gcl.variant_part_idx_type(v_key, i) == 'CDS':
						if gcl.variant_part_idx_dir(v_key, i) == t_dir:
							t_type = t_type + '-CDS'
							ds_CDS_idx = i
							ds_CDS_name = gcl.variant_part_idx_name(v_key, i)
							break
						else:
							t_type = t_type + '-rCDS'
							ds_CDS_idx = i
							ds_CDS_name = gcl.variant_part_idx_name(v_key, i)
							break
					# Move to next (previous) element
					if t_dir == 'F':
						i += 1
					else:
						i -= 1
				# Check that valid type was found
				t_details['us_CDS_idx'] = us_CDS_idx
				t_details['ds_CDS_idx'] = ds_CDS_idx
				t_details['us_Promoter_idx'] = us_Promoter_idx
				t_details['ds_Promoter_idx'] = ds_Promoter_idx
				t_details['us_Terminator_idx'] = us_Terminator_idx
				t_details['ds_Terminator_idx'] = ds_Terminator_idx
				t_details['us_CDS_name'] = us_CDS_name
				t_details['ds_CDS_name'] = ds_CDS_name
				t_details['us_Promoter_name'] = us_Promoter_name
				t_details['ds_Promoter_name'] = ds_Promoter_name
				t_details['us_Terminator_name'] = us_Terminator_name
				t_details['ds_Terminator_name'] = ds_Terminator_name
				if t_include == True:
					t_data.append([v_key, t_inst, t_name, t_type, t_details])
				else:
					t_bad_data.append([v_key, t_inst, t_name, t_type, t_details])
		return t_data, t_bad_data

##################################################################################
# TERMINATOR STRENGTH CALCULATIONS
##################################################################################

class TerminationFromTrace ():

	def __init__ (self, fwd_skip=0, fwd_len=30, rev_skip=0, rev_len=30):
		self.fwd_skip = fwd_skip
		self.fwd_len = fwd_len
		self.rev_skip = rev_skip
		self.rev_len = rev_len

	def append_strength (self, gcl, traces, t_data):
		# For each terminator calcuate the strength and append to t_data element
		for t_inst in t_data:
			v_key = t_inst[0]
			t_idx = t_inst[1]
			t_details = t_inst[4]
			trace = traces[v_key]
			t_strength = self.terminator_strength(gcl, v_key, t_idx, 
				                                  t_details['ds_CDS_idx'], 
				                                  t_details['us_CDS_idx'],
				                                  traces[v_key])
			t_inst.append(t_strength)

	def terminator_strength (self, gcl, variant, t_idx, ds_CDS_idx, us_CDS_idx, trace):
		"""Calculate an absolute terminator strength based on transcriptomic data.
		"""
		var_data = gcl.variant_data(variant)
		pre_tx = []
		post_tx = []
		# Find the direction of the TU
		step = 1
		if gcl.variant_part_idx_dir(variant, t_idx) != 'F':
			step = -1
		t_start_idx = gcl.variant_part_idx_seq_idx(variant, t_idx)
		# Strand alters where data comes from, extract relevant information
		if step == 1:
			# +ve strand
			t_end_idx = t_start_idx + gcl.variant_part_idx_seq_len(variant, t_idx)
			pre_tx = trace[0][t_start_idx-self.rev_skip-self.rev_len:t_start_idx-self.rev_skip]
			post_tx = trace[0][t_end_idx+self.fwd_skip:t_end_idx+self.fwd_skip+self.fwd_len]
		else:
			# -ve strand
			t_end_idx = t_start_idx - gcl.variant_part_idx_seq_len(variant, t_idx)
			pre_tx = trace[1][t_start_idx+self.rev_skip:t_start_idx+self.rev_skip+self.rev_len]
			post_tx = trace[1][t_end_idx-self.fwd_skip-self.fwd_len:t_end_idx-self.fwd_skip]
		# Terminator strength is the % termination (median used to smooth noise)
		t_strength = None
		if np.median(pre_tx) != 0.0:
			t_strength = 100.0-((100.0/np.median(pre_tx)) * np.median(post_tx))
		return [np.median(pre_tx), np.median(post_tx), t_strength]

class TerminationFromRSEM ():

	def __init__ (self, rsem_filename):
		rsem_data = {}
		file_reader = csv.reader(open(rsem_filename, 'rU'), delimiter=',')
		# Process header for column names (used in dict returned)
		header = next(file_reader)
		# Process each line
		for row in file_reader:
			if len(row) == 6:
				t_variant = row[0]
				t_idx = int(row[1])
				t_name = row[2]
				t_tx_in = float(row[3])
				t_tx_out = float(row[4])
				t_strength = float(row[5])
				if t_name not in rsem_data.keys():
					rsem_data[t_name] = []
				rsem_data[t_name].append([t_variant, t_idx, [t_tx_in,t_tx_out], t_strength])
		self.rsem_data = rsem_data

	def append_strength (self, t_data):
		# For each terminator calcuate the strength and append to t_data element
		for t_inst in t_data:
			v_key = t_inst[0]
			t_idx = t_inst[1]
			t_name = t_inst[2]
			t_strength = self.terminator_strength(v_key, t_idx, t_name)
			t_inst.append(t_strength)

	def terminator_strength (self, variant, t_idx, t_name):
		"""Calculate an absolute promoter strength based on RSEM data.
		"""
		found_data = [-1.0, -1.0, -1.0]
		for t_inst in self.rsem_data[t_name]:
			if t_inst[0] == variant and t_inst[1] == t_idx:
				found_data = [t_inst[2][0], t_inst[2][1], t_inst[3]]
				break
		return found_data
