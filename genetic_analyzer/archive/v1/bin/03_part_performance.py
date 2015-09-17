"""
	From an eXpress run extract the promoter and terminator performance values. Also
	process the gene level data.
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
# SUPPORTING FUNCTIONS
###############################################################################

def create_dir_if_needed (directory):
	"""Create a directory if it doesn't exist
	"""
	if not os.path.exists(directory):
		os.makedirs(directory)

def load_express_gene_results (filename):
	"""Load the eXpress gene results file into a dict of dicts.
	"""
	express_data = []
	file_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process header for column names (used in dict returned)
	header = next(file_reader)
	header_map = {}
	for idx in range(len(header)):
		header_map[idx] = header[idx]
	# Process each line
	for row in file_reader:
		if len(row) == 15:
			tx_data = {}
			tx_data[header_map[2]] = int(row[2])
			tx_data[header_map[3]] = float(row[3])
			tx_data[header_map[4]] = int(row[4])
			tx_data[header_map[5]] = int(row[5])
			tx_data[header_map[6]] = float(row[6])
			tx_data[header_map[7]] = float(row[7])
			tx_data[header_map[8]] = float(row[8])
			tx_data[header_map[9]] = float(row[9])
			tx_data[header_map[10]] = float(row[10])
			tx_data[header_map[11]] = float(row[11])
			tx_data[header_map[12]] = float(row[12])
			tx_data[header_map[13]] = row[13]
			tx_data[header_map[14]] = float(row[14])
			# Extract other attributes
			id_name = row[1].split('|')
			p_idxs = []
			p_bps = []
			tx_data['name'] = row[1]
			if 'SYNTHETIC_' in row[1]:
				id_parts = id_name[0].split('_')
				if len(id_parts) >= 4: 
					tx_data['dir'] = id_parts[-3]
					p_idxs = id_parts[-2].split('-')
					p_bps = id_parts[-1].split('-')
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
			express_data.append(tx_data)
	return express_data

###############################################################################
# EXTRACT THE PERFORMANCE ESTIMATES
###############################################################################

def run_express_extract (GCL_DATA, designs_to_process, MEASURE_TO_USE, EXPRESS_DATA_PREFIX, RESULTS_PREFIX):
	"""Extract eXpress performance data from a set of designs
	"""

	EXPRS_FILENAME = 'results_norm.xprs'

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
	for design in designs_to_process:
		# Add placeholder for data
		p_perf_data[design] = {}

		# Load the data
		design_file = EXPRESS_DATA_PREFIX+design+'/'+EXPRS_FILENAME
		design_data = load_express_gene_results(design_file)

		# Find all the promoters in the transcripts (check type of end_idx)
		p_list = []
		for trans in design_data:
			if 'SYNTHETIC_' in trans['name']:
				part_idx = trans['start_idx']
				if part_idx == None:
					if part_idx not in p_list:
							p_list.append(part_idx)
				elif gcl_data.variant_part_idx_type(design, part_idx) == 'Promoter':
						if part_idx not in p_list:
							p_list.append(part_idx)

		# Process each promoter
		for promoter_idx in p_list:
			if promoter_idx != None:
				promoter_dir = gcl_data.variant_part_idx_dir(design, promoter_idx)
				tx_readthrough = []
				tx_promoter = []
				# Go through each transcript and classify if in before or after group
				for trans in design_data:
					if 'SYNTHETIC_' in trans['name']:
						if promoter_dir == 'F':
							if (trans['start_idx'] == None or trans['start_idx'] < promoter_idx) and (trans['end_idx'] > promoter_idx or trans['end_idx'] == None):
								tx_readthrough.append(trans[MEASURE_TO_USE])
							if trans['start_idx'] == promoter_idx and (trans['end_idx'] > promoter_idx or trans['end_idx'] == None):
								tx_promoter.append(trans[MEASURE_TO_USE])
						else:
							if trans['start_idx'] == promoter_idx and (trans['end_idx'] < promoter_idx or trans['end_idx'] == None):
								tx_promoter.append(trans[MEASURE_TO_USE])
							if (trans['start_idx'] == None or trans['start_idx'] > promoter_idx) and (trans['end_idx'] < promoter_idx or trans['end_idx'] == None):
								tx_readthrough.append(trans[MEASURE_TO_USE])
				# Calculate strength
				p_perf_data[design][promoter_idx] = [np.sum(tx_readthrough), np.sum(tx_promoter)]

	# Write the performance data to file
	p_out = open(RESULTS_PREFIX+'promoter_perf.txt', 'w')
	# Write header
	p_out.write('variant,promoter_idx,promoter_name,tx_in('+MEASURE_TO_USE+'),tx_promoter('+MEASURE_TO_USE+')\n')
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
	for design in designs_to_process:
		# Add placeholder for data
		t_perf_data[design] = {}

		# Load the data
		design_file = EXPRESS_DATA_PREFIX+design+'/'+EXPRS_FILENAME
		design_data = load_express_gene_results(design_file)

		# Find all the terminators in the transcripts (check type of end_idx)
		t_list = []
		for trans in design_data:
			if 'SYNTHETIC_' in trans['name']:
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
				if 'SYNTHETIC_' in trans['name']:
					if term_dir == 'F':
						if (trans['start_idx'] == None or trans['start_idx'] < term_idx) and trans['end_idx'] == term_idx:
							tx_before.append(trans[MEASURE_TO_USE])
						if (trans['start_idx'] == None or trans['start_idx'] < term_idx) and (trans['end_idx'] > term_idx or trans['end_idx'] == None):
							tx_after.append(trans[MEASURE_TO_USE])
					else:
						if (trans['start_idx'] == None or trans['start_idx'] > term_idx) and trans['end_idx'] == term_idx:
							tx_before.append(trans[MEASURE_TO_USE])
						if (trans['start_idx'] == None or trans['start_idx'] > term_idx) and (trans['end_idx'] < term_idx or trans['end_idx'] == None):
							tx_after.append(trans[MEASURE_TO_USE])
			if np.sum(tx_before) == 0.0 and np.sum(tx_after) == 0.0:
				t_perf_data[design][term_idx] = [np.sum(tx_before), np.sum(tx_after), -1.0]
			else:
				t_perf_data[design][term_idx] = [np.sum(tx_before), np.sum(tx_after), 
				                                 np.sum(tx_before)/(np.sum(tx_before)+np.sum(tx_after))]

	# Write the performance data to file
	t_out = open(RESULTS_PREFIX+'terminator_perf.txt', 'w')
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

	# ---------------
	# GENE LEVEL DATA
	# ---------------

	# List of the extracted performances
	gene_perf_data = {}

	# For each design extract performance data
	for design in designs_to_process:
		# Add placeholder for data
		gene_perf_data[design] = {}

		# Load the data
		design_file = EXPRESS_DATA_PREFIX+design+'/'+EXPRS_FILENAME
		design_data = load_express_gene_results(design_file)

		for gene_data in design_data:
			cur_data = gene_data
			new_data = []
			new_data.append(cur_data['length'])
			new_data.append(cur_data['eff_length'])
			new_data.append(cur_data['tot_counts'])
			new_data.append(cur_data['uniq_counts'])
			new_data.append(cur_data['est_counts'])
			new_data.append(cur_data['eff_counts'])
			new_data.append(cur_data['ambig_distr_alpha'])
			new_data.append(cur_data['ambig_distr_beta'])
			new_data.append(cur_data['fpkm'])
			new_data.append(cur_data['fpkm_conf_low'])
			new_data.append(cur_data['fpkm_conf_high'])
			new_data.append(cur_data['solvable'])
			new_data.append(cur_data['tpm'])

			if 'SYNGENE_' in cur_data['name'] or 'SYNTHETIC_' in cur_data['name']:
				gene_id = cur_data['name']
				gene_perf_data[design][gene_id] = new_data
			else:
				gene_id = cur_data['name'].split('|')[0]
				gene_perf_data[design][gene_id] = new_data
		
	# Write the host genome data
	data_out = open(RESULTS_PREFIX+'gene_perf.txt', 'w')
	# Write header
	data_out.write('design,gene_id,length,eff_length,tot_counts,uniq_counts,est_counts,eff_counts,ambig_distr_alpha,ambig_distr_beta,fpkm,fpkm_conf_low,fpkm_conf_high,solvable,tpm\n')
	# Write strength data for host genes
	for design in sorted(gene_perf_data.keys()):
		for gene in sorted(gene_perf_data[design]):
			data_out.write(design + ',' + gene + ',' + ','.join([str(x) for x in gene_perf_data[design][gene]]) + '\n')
	data_out.close()

###############################################################################
# MAIN FUNCTION
###############################################################################

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="eXpress extract performance")
	parser.add_argument("-library",  dest="library",  required=True,  help="GeneClusterLibrary.txt", metavar="string")
	parser.add_argument("-designs", dest="designs", required=True, help="1,2,3", metavar="string")
	parser.add_argument("-measure",  dest="measure",  required=True,  help="fpkm", metavar="string")
	parser.add_argument("-data_prefix",  dest="data_prefix",  required=True,  help="/data/", metavar="string")
	parser.add_argument("-results_prefix",  dest="results_prefix",  required=True,  help="/results/", metavar="string")
	args = parser.parse_args()
	# Set global variables
	GCL_DATA = args.library
	MEASURE_TO_USE = args.measure
	EXPRESS_DATA_PREFIX = args.data_prefix
	RESULTS_PREFIX = args.results_prefix
	# Extract the designs to process
	designs_to_process = args.designs.split(',')
	run_express_extract(GCL_DATA, designs_to_process, MEASURE_TO_USE, EXPRESS_DATA_PREFIX, RESULTS_PREFIX)

if __name__ == "__main__":
	main()

