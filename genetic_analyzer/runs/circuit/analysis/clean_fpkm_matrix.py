#!/usr/bin/env python
"""
	Generate a clean FPKM matrix
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import csv

RESULTS_PREFIX = '../results'
OUTPUT_PREFIX = './fpkm_matrix'

tube_samples = ['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8']
flask_samples = ['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8']

def load_fpkm_matrix (filename):
	fpkms = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	header = next(data_reader)
	header_map = {}
	for idx in range(len(header)):
		if header[idx] != '':
			header_map[idx] = header[idx]
			fpkms[header[idx]] = {}
	# Process each line
	for row in data_reader:
		if len(row) > 2:
			for idx in range(1,len(header)):
				fpkms[header_map[idx]][row[0]] = float(row[idx])
	return fpkms

def make_locus_id_and_name (gene_tag):
	""" Tag format: 17_locusTag_ECDH10B_1337
	"""
	locus_id = ''
	gene_name = ''
	tag_parts = gene_tag.split('_locusTag_')
	if len(tag_parts) == 2:
		locus_id = tag_parts[1]
		gene_name = tag_parts[0]
	else:
		locus_id = 'SYNTHETIC_'+gene_tag
		gene_name = gene_tag
	return locus_id, gene_name

def save_reformatted_fpkms (fpkm_data, sample_order, filename_out):
	f_out = open(filename_out, 'w')
	header = ['locus_id', 'gene_name']+sample_order
	f_out.write('\t'.join(header)+'\n')
	# First process SYNTHETIC genes
	syn_gene_data = {}
	for s in sample_order:
		for g in sorted(fpkm_data[s].keys()):
			g_locus_id, g_name = make_locus_id_and_name(g)
			if g_locus_id.startswith('SYNTHETIC_'):
				if g not in syn_gene_data.keys():
					syn_gene_data[g] = [g_locus_id, g_name]
				syn_gene_data[g].append(str(fpkm_data[s][g]))
	# Then process ENDOGENOUS genes
	gene_data = {}
	for s in sample_order:
		for g in sorted(fpkm_data[s].keys()):
			g_locus_id, g_name = make_locus_id_and_name(g)
			if not g_locus_id.startswith('SYNTHETIC_'):
				if g not in gene_data.keys():
					gene_data[g] = [g_locus_id, g_name]
				gene_data[g].append(str(fpkm_data[s][g]))
	# Save to file
	for k in sorted(syn_gene_data.keys()):
		f_out.write('\t'.join(syn_gene_data[k])+'\n')
	for k in sorted(gene_data.keys()):
		f_out.write('\t'.join(gene_data[k])+'\n')
	f_out.close()

# Load the original data
fpkm_data = load_fpkm_matrix(RESULTS_PREFIX+'/fpkm.normed.matrix.txt')

# Reformat FPKMs
sample_order = tube_samples+flask_samples
save_reformatted_fpkms(fpkm_data, sample_order, OUTPUT_PREFIX+'/fpkm_data.txt')

# Reformat the counts too
# Load the original data
count_data = load_fpkm_matrix(RESULTS_PREFIX+'/counts.matrix.txt')

# Reformat FPKMs
sample_order = tube_samples+flask_samples
save_reformatted_fpkms(count_data, sample_order, OUTPUT_PREFIX+'/count_data.txt')


