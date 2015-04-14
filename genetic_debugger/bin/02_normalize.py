#!/usr/bin/env python
"""
	Normalize the data from an eXpress run by using DESeq to calculate
	factors that can then be applied.
"""
#	Copyright (C) 2014 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Required modules
import os
import csv
import subprocess
import argparse
import math

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

###############################################################################
# RUN THE WORKFLOW
###############################################################################

def run_express_normalizer (designs_to_process, RESULTS_PREFIX, normaliser='DESeq', inc_syn_genes=True,
	                        cnts_to_use='uniq_counts', len_to_use='length'):
	"""Run the normalizer on eXpress generated data
	"""
	gene_data = {}
	for design in designs_to_process:
		# Check file exists 
		if os.path.isfile(RESULTS_PREFIX+design+'/results.xprs') == True:
			# Load the mapped gene data
			gene_data[design] = load_express_gene_results(RESULTS_PREFIX+design+'/results.xprs')

	# Save the unique mapped reads to matrix
	read_out = open(RESULTS_PREFIX+'counts.txt', 'w')
	# Write header (designs)
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
			# We use the effective counts, not unique counts: 'uniq_counts'
			read_out.write('\t' + str(round(gene_data[design][gene][cnts_to_use],0)))
		read_out.write('\n')
	read_out.close()

	# Run normaliser on matrix to generate correction factors
	norm_filename = RESULTS_PREFIX+design+'/results_norm.xprs'
	if normaliser == 'DESeq':
		cmd_norm = 'Rscript normalize_DESeq.r ' + RESULTS_PREFIX 
		print 'express_normalizer.py RUNNING:', cmd_norm
		subprocess.call(cmd_norm, shell=True)
		norm_filename = RESULTS_PREFIX+design+'/results_norm_deseq.xprs'
	elif normaliser == 'edgeR':
		cmd_norm = 'Rscript normalize_edgeR.r ' + RESULTS_PREFIX 
		print 'express_normalizer.py RUNNING:', cmd_norm
		subprocess.call(cmd_norm, shell=True)
		norm_filename = RESULTS_PREFIX+design+'/results_norm_edger.xprs'
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

	for design in gene_data.keys():
		#if normaliser == 'edgeR':
		#	print design_lib_size[design], factors[design][1]
		for gene_name in gene_data[design].keys(): 
			if normaliser == 'DESeq':
				if gene_data[design][gene_name][len_to_use] == 0.0:
					gene_data[design][gene_name]['fpkm'] = 0.0
				else:
					gene_data[design][gene_name]['fpkm'] = (math.pow(10,9) * gene_data[design][gene_name][cnts_to_use])/float(design_lib_size[design]*gene_data[design][gene_name][len_to_use]*factors[design][0])
			elif normaliser == 'edgeR':
				#eff_lib_size = factors[design][1]*factors[design][0]
				if gene_data[design][gene_name][len_to_use] == 0.0: # or eff_lib_size == 0.0:
					gene_data[design][gene_name]['fpkm'] = 0.0
				else:
					gene_data[design][gene_name]['fpkm'] = (math.pow(10,9) * gene_data[design][gene_name][cnts_to_use])/float(design_lib_size[design]*gene_data[design][gene_name][len_to_use]*factors[design][0])
					#gene_data[design][gene_name]['fpkm'] = (math.pow(10,9) * gene_data[design][gene_name][cnts_to_use])/float(eff_lib_size*gene_data[design][gene_name][len_to_use])
	# Write the updated data to file
	for design in gene_data.keys():
		norm_filename = RESULTS_PREFIX+design+'/results_norm.xprs'
		if normaliser == 'DESeq':
			norm_filename = RESULTS_PREFIX+design+'/results_norm_deseq.xprs'
		elif normaliser == 'edgeR':
			norm_filename = RESULTS_PREFIX+design+'/results_norm_edger.xprs'

		norm_out = open(norm_filename, 'w')
		norm_head = ['bundle_id','target_id', 'length', 'eff_length', 
		             'tot_counts', 'uniq_counts', 'est_counts', 'eff_counts', 
		             'ambig_distr_alpha', 'ambig_distr_beta', 'fpkm', 
		             'fpkm_conf_low', 'fpkm_conf_high', 'solvable', 'tpm']
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



# Some tests
designs_to_process = ['1_1', '2_1', '3_1', '4_1', '5_1', '6_1', '7_1', '8_1', '1_2', '2_2', '3_2', '4_2', '5_2', '6_2', '7_2', '8_2']
run_express_normalizer(designs_to_process, '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/circuit/exp_gene/',
	                   normaliser='DESeq', inc_syn_genes=True, cnts_to_use='uniq_counts', len_to_use='length')
run_express_normalizer(designs_to_process, '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/circuit/exp_transcript/',
	                   normaliser='DESeq', inc_syn_genes=False, cnts_to_use='eff_counts', len_to_use='eff_length')
run_express_normalizer(designs_to_process, '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/circuit/exp_gene/',
	                   normaliser='edgeR', inc_syn_genes=True, cnts_to_use='uniq_counts', len_to_use='length')
run_express_normalizer(designs_to_process, '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/circuit/exp_transcript/',
	                   normaliser='edgeR', inc_syn_genes=False, cnts_to_use='eff_counts', len_to_use='eff_length')

# Which designs to create estimates for (85 to process all designs)
designs_to_process = [str(x) for x in range(1,85)]
# Remove designs with sequencing errors
temp_designs = []
for idx in range(len(designs_to_process)):
	if designs_to_process[idx] not in ['17', '18', '33', '45', '57', '69', '76', '81']:
		temp_designs.append(designs_to_process[idx])
designs_to_process = temp_designs
run_express_normalizer(designs_to_process,
	                   '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/stata/exp_gene/',
	                   normaliser='DESeq', inc_syn_genes=True, cnts_to_use='uniq_counts', len_to_use='length')
run_express_normalizer(designs_to_process,
	                   '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/stata/exp_transcript/',
	                   normaliser='DESeq', inc_syn_genes=False, cnts_to_use='eff_counts', len_to_use='eff_length')
# cnts_to_use='eff_counts', len_to_use='eff_length'
run_express_normalizer(designs_to_process,
	                   '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/stata/exp_gene/',
	                   normaliser='edgeR', inc_syn_genes=True, cnts_to_use='uniq_counts', len_to_use='length')
run_express_normalizer(designs_to_process,
	                   '/Users/Tom/Dropbox/Research/projects/MIT_Voigt_Lab/rnaseq_workflow/results/stata/exp_transcript/',
	                   normaliser='edgeR', inc_syn_genes=False, cnts_to_use='eff_counts', len_to_use='eff_length')

###############################################################################
# MAIN FUNCTION
###############################################################################

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="eXpress run")
	parser.add_argument("-designs",  dest="designs",  required=True,  help="1,2,3", metavar="string")
	parser.add_argument("-fastq_prefix",  dest="fastq_prefix",  required=True,  help="/fastq/", metavar="string")
	parser.add_argument("-data_prefix",  dest="data_prefix",  required=True,  help="/data/", metavar="string")
	parser.add_argument("-results_gene_prefix",  dest="results_gene_prefix",  required=True,  help="/results_gene/", metavar="string")
	parser.add_argument("-results_prefix",  dest="results_prefix",  required=True,  help="/results/", metavar="string")
	parser.add_argument("-normaliser",  dest="normaliser",  required=False,  help="edgeR", metavar="string")
	args = parser.parse_args()
	# Set global variables
	FASTQ_PREFIX = args.fastq_prefix
	DATA_PREFIX = args.data_prefix
	RESULTS_GENE_PREFIX = args.results_gene_prefix
	RESULTS_PREFIX = args.results_prefix
	NORMER = 'edgeR'
	if args.normaliser:
		NORMER = args.normaliser
	# Extract the designs to process and run
	designs_to_process = args.designs.split(',')
	# Run RSEM commands
	run_express_normalizer(designs_to_process, RESULTS_GENE_PREFIX, RESULTS_PREFIX, normaliser=NORMER)

#if __name__ == "__main__":
#	main()
