#!/usr/bin/env python
"""
	Generate a clean GFF file
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import csv

def load_fasta (filename):
	fasta_data = {}
	f_in = open(filename, 'rU')
	line_data = f_in.readlines()
	f_in.close()
	cur_chrom = None
	cur_seq = ''
	for line in line_data:
		cur_line = line.lstrip().rstrip()
		if len(cur_line) > 0:
			if cur_line[0] == '>':
				# End prev seq
				if cur_chrom != None:
					fasta_data[cur_chrom] = cur_seq
				chrom_parts = cur_line[1:].split('|')
				cur_chrom = chrom_parts[3]
				cur_seq = ''
			else:
				cur_seq += cur_line
	if cur_chrom != None:
		fasta_data[cur_chrom] = cur_seq
	return fasta_data

def save_reformatted_fasta (fasta_data, filename_out):
	f_out = open(filename_out, 'w')
	# Save to file
	for chrom in sorted(fasta_data.keys()):
		f_out.write('>'+chrom+'\n')
		f_out.write(fasta_data[chrom]+'\n\n')
	f_out.close()

original_fasta = load_fasta('Koxytoca5al_org.fasta')
save_reformatted_fasta(original_fasta, 'Koxytoca5al.fasta')
