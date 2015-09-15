#!/usr/bin/env python
"""
	Normalize the data from an eXpress.
"""
#	Copyright (C) 2014 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Required modules
import os
import csv
import argparse
import math
import gene_cluster_library as gcl

###############################################################################
# SUPPORTING FUNCTIONS
###############################################################################

def load_express_gene_results (filename):
	"""Load the eXpress gene results file into a dict of dicts.
	"""
	express_data = {}
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
			tx_data[header_map[0]] = int(row[0])
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
			# Add the data element
			express_data[row[1]] = tx_data
	return express_data

def load_factors (filename):
	factor_data = {}
	file_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process header for column names (used in dict returned)
	header = next(file_reader)
	for row in file_reader:
		if len(row) >= 2:
			des = row[0]
			if row[0][0] == 'X':
				des = row[0][1:]
			row_data = [float(x) for x in row[1:]]
			factor_data[des] = row_data
	return factor_data

def cds_names_in_idx_range (gcl_data, design, idx_start, idx_end, gene_dir):
	syngenes = {} # Gene name as key, length as value
	# Check direction
	if gene_dir != 'F':
		tmp = idx_start
		idx_start = idx_end
		idx_end = tmp
	# Make sure indexes are correct
	if idx_start == None:
		idx_start = 0
	if idx_end == None:
		idx_end = len(gcl_data.variants[design]['part_list'])
	# Output all CDSs
	for cur_idx in range(idx_start, idx_end+1):
		cur_part_name = gcl_data.variants[design]['part_list'][cur_idx]['part_name']
		if ( gcl_data.variant_part_idx_type(design, cur_idx) == 'CDS' and 
			 gcl_data.variant_part_idx_dir(design, cur_idx) == gene_dir ):
			new_name = 'SYNGENE_'+cur_part_name+'_'+cur_idx
			syngenes[new_name] = gcl_data.variant_part_idx_seq_len(design, cur_idx)
	return syngenes

###############################################################################
# RUN THE WORKFLOW
###############################################################################

def run_normalizer (GCL_DATA, designs_to_process, RESULTS_PREFIX, normaliser='edgeR', inc_syn_genes=True,
	                cnts_to_use='est_counts', len_to_use='length'):
	"""Run the normalizer (cnts_to_use can be: est_counts, eff_counts, uniq_counts)
	"""

	# All header items
	norm_head = ['bundle_id','target_id', 'length', 'eff_length', 
		         'tot_counts', 'uniq_counts', 'est_counts', 'eff_counts', 
		         'ambig_distr_alpha', 'ambig_distr_beta', 'fpkm', 
		         'fpkm_conf_low', 'fpkm_conf_high', 'solvable', 'tpm']

	gene_data = {}
	for design in designs_to_process:
		# Check file exists 
		if os.path.isfile(RESULTS_PREFIX+design+'/results.xprs') == True:
			# Load the mapped gene data
			gene_data[design] = load_express_gene_results(RESULTS_PREFIX+design+'/results.xprs')

	# Save the unique mapped reads to matrix
	read_out = open(RESULTS_PREFIX+'counts.txt', 'w')
	# Write header (designs)
	read_out.write('transcript')
	designs_found = sorted(gene_data.keys())
	for design in designs_found:
		read_out.write('\t' + design)
	read_out.write('\n')
	gene_ids = gene_data[designs_found[0]].keys()
	if inc_syn_genes == False:
		gene_ids = [x for x in gene_ids if 'SYNTHETIC_' not in x]
	for gene in sorted(gene_ids):
		read_out.write(gene)
		for design in designs_found:
			read_out.write('\t' + str(round(gene_data[design][gene][cnts_to_use],0)))
		read_out.write('\n')
	read_out.close()

	# Run normaliser on matrix to generate correction factors
	if normaliser == 'DESeq':
		cmd_norm = 'Rscript normalize_DESeq.r ' + RESULTS_PREFIX 
		print 'express_normalizer.py RUNNING:', cmd_norm
		subprocess.call(cmd_norm, shell=True)
	elif normaliser == 'edgeR':
		cmd_norm = 'Rscript normalize_edgeR.r ' + RESULTS_PREFIX 
		print 'express_normalizer.py RUNNING:', cmd_norm
		subprocess.call(cmd_norm, shell=True)
	else:
		print 'express_normalizer.py ERROR: normalizer type not known, ', normaliser
		return

	# Apply factors to fpkm values from eXpress and save normalised output (results_norm.xprs)
	factors = {}
	if normaliser == 'DESeq':
		factors = load_factors(RESULTS_PREFIX+'correction_factors_deseq.txt')
	elif normaliser == 'edgeR':
		factors = load_factors(RESULTS_PREFIX+'correction_factors_edger.txt')

	design_lib_size = {}
	for design in gene_data.keys():
		design_lib_size[design] = 0.0
		for gene_name in gene_data[design].keys():
			design_lib_size[design] += gene_data[design][gene_name][cnts_to_use]

	# Generate counts for synthetic genes aswell

	# Load the design data (needed to link genes to part idxs)
	gcl_data = gcl.GeneClusterLibrary()
	gcl_data.load(GCL_DATA)

	# Cycle through each SYNTHETIC transcript and see if present
	for design in gene_data.keys():
		for gene_name in gene_data[design].keys():
			if 'SYNTHETIC_' in gene_name:
				# Extract start and end indexes
				idx_start = 0
				idx_end = 0
				id_parts = gene_name.split('_')
				gene_dir = id_parts[-3]
				p_idxs = id_parts[-2].split('-')
				if p_idxs[0] == 'None':
					idx_start = None
				else:
					idx_start = int(p_idxs[0])
				if p_idxs[1] == 'None':
					idx_end = None
				else:
					idx_end = int(p_idxs[1])
				syngenes = cds_names_in_idx_range(gcl_data, design, idx_start, idx_end, gene_dir)
				# Update the counts of each gene
				for syngene in syngenes.keys():
					if syngene not in gene_data[design].keys():
						gene_data[design][syngene] = {}
						for h in norm_head:
							if h == 'target_id':
								gene_data[design][syngene][h] = syngene
							elif h == 'solvable':
								gene_data[design][syngene][h] = ''
							elif h == 'length' or h == 'eff_length':
								gene_data[design][syngene][h] = syngenes[syngene]
							else:
								gene_data[design][syngene][h] = 0
					# Add counts
					gene_data[design][syngene]['tot_counts'] += gene_data[design][gene_name]['tot_counts']
					gene_data[design][syngene]['uniq_counts'] += gene_data[design][gene_name]['uniq_counts']
					gene_data[design][syngene]['est_counts'] += gene_data[design][gene_name]['est_counts']
					gene_data[design][syngene]['eff_counts'] += gene_data[design][gene_name]['eff_counts']

	for design in gene_data.keys():
		for gene_name in gene_data[design].keys(): 
			if normaliser == 'DESeq':
				if gene_data[design][gene_name][len_to_use] == 0.0:
					gene_data[design][gene_name]['fpkm'] = 0.0
				else:
					gene_data[design][gene_name]['fpkm'] = (math.pow(10,9) * gene_data[design][gene_name][cnts_to_use])/float(design_lib_size[design]*gene_data[design][gene_name][len_to_use]*factors[design][0])
			elif normaliser == 'edgeR':
				if gene_data[design][gene_name][len_to_use] == 0.0:
					gene_data[design][gene_name]['fpkm'] = 0.0
				else:
					gene_data[design][gene_name]['fpkm'] = (math.pow(10,9) * gene_data[design][gene_name][cnts_to_use])/float(design_lib_size[design]*gene_data[design][gene_name][len_to_use]*factors[design][0])
					
	# Write the updated data to file
	for design in gene_data.keys():
		norm_filename = RESULTS_PREFIX+design+'/results_norm.xprs'
		norm_out = open(norm_filename, 'w')
		norm_out.write('\t'.join(norm_head) + '\n')
		for gene_name in sorted(gene_data[design].keys()):
			first_f = True
			for f in norm_head:
				if first_f == False:
					norm_out.write('\t')
				else:
					first_f = False
				if f == 'target_id':
					norm_out.write(gene_name)
				else:
					norm_out.write(str(gene_data[design][gene_name][f]))
			norm_out.write('\n')
		norm_out.close()

###############################################################################
# MAIN FUNCTION
###############################################################################

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="normalize")
	parser.add_argument("-library",  dest="library",  required=True,  help="GeneClusterLibrary.txt", metavar="string")
	parser.add_argument("-designs", dest="designs", required=True, help="1,2,3", metavar="string")
	parser.add_argument("-results_prefix", dest="results_prefix", required=True, help="/results/", metavar="string")
	parser.add_argument("-normaliser", dest="normaliser", required=False, help="edgeR", metavar="string")
	args = parser.parse_args()
	# Set global variables
	GCL_DATA = args.library
	RESULTS_PREFIX = args.results_prefix
	valid_normalisers = ['edgeR', 'DESeq']
	NORMER = 'edgeR'
	if args.normaliser:
		NORMER = args.normaliser
	# Extract the designs to process
	designs_to_process = args.designs.split(',')
	run_express_normalizer(GCL_DATA, designs_to_process, RESULTS_PREFIX, normaliser=NORMER)

if __name__ == "__main__":
	main()
