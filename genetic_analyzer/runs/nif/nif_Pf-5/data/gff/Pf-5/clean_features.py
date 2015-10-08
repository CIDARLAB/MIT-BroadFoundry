#!/usr/bin/env python
"""
	Generate a clean GFF file
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import csv

def load_feature_data (filename):
	feature_data = []
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process each line
	start_bp = 0
	end_bp = 0
	gene_name = 'None'
	locus_tag = 'None'
	gene_dir = '+'
	for row in data_reader:
		if len(row) >= 3:
			if row[0] != '' and row[1] != '' and row[2] == 'gene':
				# We have a new gene, add old one
				if locus_tag != 'None':
					tag = gene_name.replace(' ', '')+'_locusTag_'+locus_tag
					feature_data.append([tag, start_bp, end_bp, gene_dir])
				# Update next features position
				start_bp = int(row[0])
				end_bp = int(row[1])
				if start_bp > end_bp:
					temp_end_bp = end_bp
					end_bp = start_bp
					start_bp = temp_end_bp
					gene_dir = '-'
				gene_name = 'None'
				locus_tag = 'None'
			else:
				if len(row) >= 4:
					if row[3] == 'gene':
						gene_name = row[4]
					if row[3] == 'locus_tag':
						locus_tag = row[4]
	if locus_tag != 'None':
		tag = gene_name.replace(' ', '')+'_locusTag_'+locus_tag
		feature_data.append([tag, start_bp, end_bp, gene_dir])
	return feature_data

def save_reformatted_deg_data (feature_data, chrom, filename_out):
	f_out = open(filename_out, 'w')
	# Save to file
	for el in feature_data:
		output_line = [chrom, 'RefSeq', 'gene', str(el[1]), str(el[2]), '.', el[3], '.', 'Name='+el[0]]
		f_out.write('\t'.join(output_line)+'\n')
	f_out.close()

feature_data = load_feature_data('./genbank_features.txt')
save_reformatted_deg_data(feature_data, 'Pseudomonas_pf5', 'Pseudomonas_pf5.gff')
save_reformatted_deg_data(feature_data, 'Pseudomonas_pf5_weakT7', 'Pseudomonas_pf5_weakT7.gff')
save_reformatted_deg_data(feature_data, 'Pseudomonas_pf5_stongT7', 'Pseudomonas_pf5_stongT7.gff')
