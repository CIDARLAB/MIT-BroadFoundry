#!/usr/bin/env python
"""
Reformats the GO association files for use with goatools.
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import csv

def load_go_data (go_filename, gene_idx, go_term_idx):
	go_data = {}
	# Open file and ignore headers
	go_reader = csv.reader(open(go_filename, 'rb'), delimiter='\t')
	# Load all the data
	for row in go_reader:
		if len(row) > 1:
			if row[gene_idx] not in go_data.keys():
				go_data[row[gene_idx]] = []
			go_data[row[gene_idx]].append(row[go_term_idx])
	return go_data

def write_goatool_format (go_data, filename_out):
	f_out = open(filename_out, 'w')
	for g in go_data.keys():
		f_out.write(g + '\t' + ';'.join(go_data[g]) + '\n')
	f_out.close()

def write_population_gene_list (go_data, filename_out):
	f_out = open(filename_out, 'w')
	for g in go_data.keys():
		f_out.write(g + '\n')
	f_out.close()

go_data = load_go_data('gene_association.ecocyc', 2, 4)
write_goatool_format(go_data, 'ecoli_go_association.txt')
write_population_gene_list(go_data, 'ecoli_go_population.txt')
