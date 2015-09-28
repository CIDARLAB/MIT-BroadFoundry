#!/usr/bin/env python
"""
	Generate a clean DEG data files
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import csv

RESULTS_PREFIX = '../results'
OUTPUT_PREFIX = './deg_analysis'

def load_deg_data (filename):
	deg_data = []
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) > 2:
			locus_id, gene_name = make_locus_id_and_name(row[0])
			deg_data.append([locus_id, gene_name]+row[1:])
	return deg_data

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

def save_reformatted_deg_data (deg_data, filename_out):
	f_out = open(filename_out, 'w')
	header = ['locus_id', 'gene_name', 'log2(FC)', 'log2(CPM)', 'p_value', 'fdr']
	f_out.write('\t'.join(header)+'\n')
	# Save to file
	for el in deg_data:
		f_out.write('\t'.join(el)+'\n')
	f_out.close()

# Reformat all studies
deg_data = load_deg_data(RESULTS_PREFIX+'/broken_flask_vs_tube.de.analysis.txt')
save_reformatted_deg_data(deg_data, OUTPUT_PREFIX+'/broken_flask_vs_tube_de_data.txt')

deg_data = load_deg_data(RESULTS_PREFIX+'/flask_vs_tube.de.analysis.txt')
save_reformatted_deg_data(deg_data, OUTPUT_PREFIX+'/flask_vs_tube_de_data.txt')

deg_data = load_deg_data(RESULTS_PREFIX+'/ara_comp_tube.de.analysis.txt')
save_reformatted_deg_data(deg_data, OUTPUT_PREFIX+'/ara_comp_tube_de_data.txt')

deg_data = load_deg_data(RESULTS_PREFIX+'/atc_comp_tube.de.analysis.txt')
save_reformatted_deg_data(deg_data, OUTPUT_PREFIX+'/atc_comp_tube_de_data.txt')

deg_data = load_deg_data(RESULTS_PREFIX+'/iptg_comp_tube.de.analysis.txt')
save_reformatted_deg_data(deg_data, OUTPUT_PREFIX+'/iptg_comp_tube_de_data.txt')
