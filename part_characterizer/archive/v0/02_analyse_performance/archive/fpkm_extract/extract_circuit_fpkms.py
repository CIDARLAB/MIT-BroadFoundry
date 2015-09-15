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
			#tx_data['name'] = row[1]
			if 'SYNTHETIC_' in row[1]:
				id_parts = id_name[0].split('_')
				if len(id_parts) == 5:
					tx_data['dir'] = id_parts[2]
					p_idxs = id_parts[3].split('-')
					p_bps = id_parts[4].split('-')
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
				express_data[row[1]] = tx_data['fpkm']
			else:
				express_data[row[1]] = tx_data['fpkm']
	return express_data

###############################################################################
# EXTRACT THE PERFORMANCE ESTIMATES
###############################################################################

###############################################################################
# PARAMETERS
###############################################################################

# Inputs:
# 1. GeneClusterLibrary of design
# 2. Location to eXpress data (both transcript and gene level)
# 3. Location to save results to

# Which data to process?


# CIRCUIT WORKING/BROKEN
LIBRARY_TO_PROCESS = 'circuit'
#GCL_DATA = '../data/'+LIBRARY_TO_PROCESS+'/circuit_library.txt'

# Measure to use to calculate performance
MEASURE_TO_USE = 'fpkm' # tpm/fpkm

# Standard locations
EXPRESS_DATA_PREFIX = '../results/'+LIBRARY_TO_PROCESS+'/exp_transcript/'
EXPRESS_GENE_DATA_PREFIX = '../results/'+LIBRARY_TO_PROCESS+'/exp_gene/'
RESULTS_PREFIX = '../results/'+LIBRARY_TO_PROCESS+'/part_perf/'

file_to_use = 'results_norm_edger.xprs'

data_1 = {}
data_1['1'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'1_1/'+file_to_use)
data_1['2'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'2_1/'+file_to_use)
data_1['3'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'3_1/'+file_to_use)
data_1['4'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'4_1/'+file_to_use)
data_1['5'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'5_1/'+file_to_use)
data_1['6'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'6_1/'+file_to_use)
data_1['7'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'7_1/'+file_to_use)
data_1['8'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'8_1/'+file_to_use)

data_2 = {}
data_2['1'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'1_2/'+file_to_use)
data_2['2'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'2_2/'+file_to_use)
data_2['3'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'3_2/'+file_to_use)
data_2['4'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'4_2/'+file_to_use)
data_2['5'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'5_2/'+file_to_use)
data_2['6'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'6_2/'+file_to_use)
data_2['7'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'7_2/'+file_to_use)
data_2['8'] = load_express_gene_results(EXPRESS_GENE_DATA_PREFIX+'8_2/'+file_to_use)


#######################################################################
# CIRCUIT DATA
#######################################################################

cir_gene_list = ['AmtR', 'LitR', 'BM3R1', 'PhlF', 'SrpR', 'YFP']

out_data = {}

for i in [str(x) for x in range(1,9)]:
	d1 = data_1[i]
	d2 = data_2[i]
	out_data[i] = []
	for gene_id in d1.keys():
		for g in cir_gene_list:
			if g in gene_id:
				# Process
				out_data[i].append([g, d1[gene_id], d2[gene_id]])

# Save resuts to file for easy inclusion in further analysis
rf_data = {}
for g in cir_gene_list:
	rf_data[g] = {}
	for d in out_data.keys():
		for i in range(len(out_data[d])):
			if out_data[d][i][0] == g:
				rf_data[g][d]=out_data[d][i][1:]

f_out = open('circuit_fpkm_data.txt', 'w')
for gene in sorted(rf_data.keys()):
	f_out.write(gene + '\t')
	for design in sorted(rf_data[gene].keys()):
		f_out.write(str(rf_data[gene][design][0]) + '\t')
	for design in sorted(rf_data[gene].keys()):
		f_out.write(str(rf_data[gene][design][1]) + '\t')
	f_out.write('\n')
f_out.close()


#######################################################################
# GENOME DATA
#######################################################################



gene_list = []

out_data = {}

for i in [str(x) for x in range(1,9)]:
	d1 = data_1[i]
	d2 = data_2[i]
	out_data[i] = []
	for gene_id in d1.keys():
		cir_gene = False
		for g in cir_gene_list:
			if g in gene_id:
				cir_gene = True
		if cir_gene == False:
			gene_parts = gene_id.split('|')
			gene_list.append(gene_parts[0]+'\t'+gene_parts[1])
			out_data[i].append([gene_parts[0]+'\t'+gene_parts[1], d1[gene_id], d2[gene_id]])

# Save resuts to file for easy inclusion in further analysis
rf_data = {}
for g in gene_list:
	rf_data[g] = {}
	for d in out_data.keys():
		for i in range(len(out_data[d])):
			if out_data[d][i][0] == g:
				rf_data[g][d]=out_data[d][i][1:]

f_out = open('genome_fpkm_data.txt', 'w')
for gene in sorted(rf_data.keys()):
	f_out.write(gene + '\t')
	for design in sorted(rf_data[gene].keys()):
		f_out.write(str(rf_data[gene][design][0]) + '\t')
	for design in sorted(rf_data[gene].keys()):
		f_out.write(str(rf_data[gene][design][1]) + '\t')
	f_out.write('\n')
f_out.close()










