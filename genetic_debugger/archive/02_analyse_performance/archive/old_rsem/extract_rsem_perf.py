"""
	From an RSEM run extract the promoter and terminator performance values
"""
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Required modules
import os
import csv
import numpy as np
import gene_cluster_library as gcl

###############################################################################
# PARAMETERS
###############################################################################

# Which data to process?

# STATA LIBRARY
#LIBRARY_TO_PROCESS = 'stata_trunc'
#GCL_DATA = '../data/'+LIBRARY_TO_PROCESS+'/clean_nif_stata_library.txt'

# CIRCUIT WORKING/BROKEN
LIBRARY_TO_PROCESS = 'circuit_trunc'
GCL_DATA = '../data/'+LIBRARY_TO_PROCESS+'/full_circuit_lib.txt'

# Measure to use to calculate performance
MEASURE_TO_USE = 'FPKM'

# Standard locations
RSEM_DATA_PREFIX = '../results/rsem_'+LIBRARY_TO_PROCESS+'/rsem_exps/'
RESULTS_PREFIX = '../results/rsem_'+LIBRARY_TO_PROCESS+'/rsem_perf/'

###############################################################################
# SUPPORTING FUNCTIONS
###############################################################################

def create_dir_if_needed (directory):
	"""Create a directory if it doesn't exist
	"""
	if not os.path.exists(directory):
		os.makedirs(directory)

def load_rsem_gene_results (filename):
	"""Load the RSEM gene results file into a dict of dicts.
	"""
	rsem_data = []
	file_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process header for column names (used in dict returned)
	header = next(file_reader)
	header_map = {}
	for idx in range(len(header)):
		header_map[idx] = header[idx]
	# Process each line
	for row in file_reader:
		if len(row) == 7:
			tx_data = {}
			tx_data[header_map[2]] = float(row[2])
			tx_data[header_map[3]] = float(row[3])
			tx_data[header_map[4]] = float(row[4])
			tx_data[header_map[5]] = float(row[5])
			tx_data[header_map[6]] = float(row[6])
			# Extract other attributes
			id_parts = row[0].split('_')
			tx_data['dir'] = id_parts[1]
			p_idxs = id_parts[2].split('-')
			p_bps = id_parts[3].split('-')
			if p_idxs[0] == 'None':
				tx_data['start_idx'] = None
				tx_data['start_bp'] = None
			else:
				tx_data['start_idx'] = int(p_idxs[0])
				tx_data['start_bp'] = int(p_bps[0])
			if p_idxs[1] == 'None':
				tx_data['end_idx'] = None
				tx_data['end_bp'] = None
			else:
				tx_data['end_idx'] = int(p_idxs[1])
				tx_data['end_bp'] = int(p_bps[1])
		# Add the data element
		rsem_data.append(tx_data)
	return rsem_data

###############################################################################
# EXTRACT THE PERFORMANCE ESTIMATES
###############################################################################

# Make sure we have somewhere to save our results
create_dir_if_needed(RESULTS_PREFIX)

# Load the design data (needed to link promoter/terminator names to part idxs)
gcl_data = gcl.GeneClusterLibrary()
gcl_data.load(GCL_DATA)

# -----------------
# PROMOTER ANALYSIS
# -----------------

# List of the extracted performances
p_perf_data = {}

# For each design extract performance data
for design in gcl_data.variants.keys():
	# Add placeholder for data
	p_perf_data[design] = {}

	# Load the data
	design_file = RSEM_DATA_PREFIX+'rsem_exp_'+design+'.genes.results'
	design_data = load_rsem_gene_results(design_file)

	# Find all the promoters in the transcripts (check type of end_idx)
	p_list = []
	for trans in design_data:
		part_idx = trans['start_idx']
		if part_idx != None and gcl_data.variant_part_idx_type(design, part_idx) == 'Promoter':
				if part_idx not in p_list:
					p_list.append(part_idx)

	# Process each promoter
	for promoter_idx in p_list:
		promoter_dir = gcl_data.variant_part_idx_dir(design, promoter_idx)
		tx_readthrough = []
		tx_promoter = []
		# Go through each transcript and classify if in before or after group
		for trans in design_data:
			if promoter_dir == 'F':
				if trans['start_idx'] < promoter_idx and trans['end_idx'] >= promoter_idx:
					tx_readthrough.append(trans[MEASURE_TO_USE])
				if trans['start_idx'] == promoter_idx and trans['end_idx'] > promoter_idx:
					tx_promoter.append(trans[MEASURE_TO_USE])
			else:
				if trans['start_idx'] == promoter_idx and trans['end_idx'] < promoter_idx:
					tx_promoter.append(trans[MEASURE_TO_USE])
				if trans['start_idx'] > promoter_idx and trans['end_idx'] < promoter_idx:
					tx_readthrough.append(trans[MEASURE_TO_USE])
		# Calculate strength
		p_perf_data[design][promoter_idx] = [np.sum(tx_readthrough), np.sum(tx_promoter)]

# Write the performance data to file
p_out = open(RESULTS_PREFIX+'promoter_rsem_pref.txt', 'w')
# Write header
p_out.write('variant,promoter_idx,promoter_name,tx_readthrough('+MEASURE_TO_USE+'),tx_promoter('+MEASURE_TO_USE+')\n')
# Write strength data
for design in sorted(p_perf_data.keys()):
	for promoter_idx in p_perf_data[design]:
		data_out = [design, str(promoter_idx), gcl_data.variant_part_idx_name(design, promoter_idx),
		            str(p_perf_data[design][promoter_idx][0]), str(p_perf_data[design][promoter_idx][1])+'\n']
		p_out.write(','.join(data_out))
p_out.close()

# -------------------
# TERMINATOR ANALYSIS
# -------------------

# List of the extracted performances
t_perf_data = {}

# For each design extract performance data
for design in gcl_data.variants.keys():
	# Add placeholder for data
	t_perf_data[design] = {}

	# Load the data
	design_file = RSEM_DATA_PREFIX+'rsem_exp_'+design+'.genes.results'
	design_data = load_rsem_gene_results(design_file)

	# Find all the terminators in the transcripts (check type of end_idx)
	t_list = []
	for trans in design_data:
		part_idx = trans['end_idx']
		if part_idx != None and gcl_data.variant_part_idx_type(design, part_idx) == 'Terminator':
				if part_idx not in t_list:
					t_list.append(part_idx)

	# Process each terminator
	for term_idx in t_list:
		term_dir = gcl_data.variant_part_idx_dir(design, term_idx)
		tx_before = []
		tx_after = []
		# Go through each transcript and classify if in before or after group
		for trans in design_data:
			if term_dir == 'F':
				if trans['start_idx'] < term_idx and trans['end_idx'] == term_idx:
					tx_before.append(trans[MEASURE_TO_USE])
				if trans['start_idx'] < term_idx and trans['end_idx'] > term_idx:
					tx_after.append(trans[MEASURE_TO_USE])
			else:
				if trans['start_idx'] > term_idx and trans['end_idx'] == term_idx:
					tx_before.append(trans[MEASURE_TO_USE])
				if trans['start_idx'] > term_idx and trans['end_idx'] < term_idx:
					tx_after.append(trans[MEASURE_TO_USE])
		if np.sum(tx_before) == 0.0 and np.sum(tx_after) == 0.0:
			t_perf_data[design][term_idx] = [np.sum(tx_before), np.sum(tx_after), -1.0]
		else:
			t_perf_data[design][term_idx] = [np.sum(tx_before), np.sum(tx_after), 
			                                 np.sum(tx_before)/(np.sum(tx_before)+np.sum(tx_after))]

# Write the performance data to file
t_out = open(RESULTS_PREFIX+'terminator_rsem_pref.txt', 'w')
# Write header
t_out.write('variant,terminator_idx,terminator_name,tx_in('+MEASURE_TO_USE+'),tx_out('+MEASURE_TO_USE+'),strength('+MEASURE_TO_USE+')\n')
# Write strength data
for design in sorted(t_perf_data.keys()):
	for terminator_idx in t_perf_data[design]:
		data_out = [design, str(terminator_idx), gcl_data.variant_part_idx_name(design, terminator_idx),
		            str(t_perf_data[design][terminator_idx][0]), str(t_perf_data[design][terminator_idx][1]), 
		            str(t_perf_data[design][terminator_idx][2])+'\n']
		t_out.write(','.join(data_out))
t_out.close()
