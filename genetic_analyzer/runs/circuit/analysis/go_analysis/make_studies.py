#!/usr/bin/env python
"""
Make the study files for GO enrichment analysis
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import csv

DEG_PREFIX = '../../results/'
OUT_PREFIX = './studies/'

def load_de_genes (filename_in, p_val=0.01, fdr=0.05):
	de_genes = {}
	with open(filename_in, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		# Skip header
		next(reader)
		for row in reader:
			if len(row) == 5:
				gene_parts = row[0].split('_locusTag_')
				if len(gene_parts) == 2:
					cur_gene = (row[0].split('_locusTag_'))[0]
					cur_fc = float(row[1])
					cur_p_val = float(row[3])
					cur_fdr = float(row[4])
					if cur_p_val <= p_val and cur_fdr <= fdr:
						de_genes[cur_gene] = [cur_fc, cur_p_val, cur_fdr]
	return de_genes

def save_study (de_genes, filename_up, filename_down):
	gene_list_up = []
	gene_list_down = []
	for g in de_genes.keys():
		if de_genes[g][0] < 0.0:
			gene_list_down.append(g)
		else:
			gene_list_up.append(g)
	f_out = open(filename_up, 'w')
	for g in gene_list_up:
		f_out.write(g+'\n')
	f_out.close()
	f_out = open(filename_down, 'w')
	for g in gene_list_down:
		f_out.write(g+'\n')
	f_out.close()

de_genes_flask_vs_tube_all = load_de_genes(DEG_PREFIX+'flask_vs_tube.de.analysis.txt', p_val=0.01, fdr=1.0)
de_genes_flask_vs_tube_678 = load_de_genes(DEG_PREFIX+'broken_flask_vs_tube.de.analysis.txt', p_val=0.01, fdr=1.0)

save_study(de_genes_flask_vs_tube_all, OUT_PREFIX+'flask_vs_tube.up.txt', OUT_PREFIX+'flask_vs_tube.down.txt')
save_study(de_genes_flask_vs_tube_678, OUT_PREFIX+'broken_flask_vs_tube.up.txt', OUT_PREFIX+'broken_flask_vs_tube.down.txt')


