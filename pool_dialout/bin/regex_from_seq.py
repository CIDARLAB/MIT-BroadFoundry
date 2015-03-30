#!/usr/bin/env python
"""
Generate reference regular expressions
======================================
	From the PCR product generate reference regular expression that has
	a fixed length in both directions and replaces Ns with barcode
	regions to extract.
"""
from __future__ import print_function, division
import os
import sys
import re
import string

__author__  = 'Thomas E. Gorochowski, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

## HELPERS ====================================================================

# Reverse complement
def revcomp(seq, trans=string.maketrans("ACGTN", "TGCAN")):
	return "".join(reversed(seq.translate(trans)))

# Replace Ns with grouped regex that can be used to extract barcode
def replace_N(seq, trim_end=True):
	new_seq = ""
	idx = 0
	while idx < len(seq):
		if seq[idx] in "ACGT":
			new_seq += seq[idx]
		else:
			n_num = 1
			while idx+1 < len(seq) and seq[idx+1] == "N":
				n_num += 1
				idx += 1
			if trim_end == False or (trim_end == True and idx < len(seq)-1):
				new_seq += "([ACGT]{" + str(n_num) + "})"
		idx += 1
	return new_seq

## MAIN =======================================================================

# Parse command line parameters
if len(sys.argv) != 4:
	print("Usage: python {} <seq fasta> <ref length> <output filename>".format(sys.argv[0]), file=sys.stderr)
	sys.exit()
seq_filename, ref_length, out_filename = sys.argv[1:]
ref_length = int(ref_length)

# Load regular expressions uniquely defining each design and compile
seqs = {}
with open(seq_filename, "rU") as seq_file:
	for header in seq_file:
		cur_seq_name = header[1:].strip()
		cur_seq = seq_file.next().strip().upper()
		seqs[cur_seq_name] = cur_seq

regex_seqs = {}
for seq_key in seqs.keys():
	cur_seq = seqs[seq_key]

	# Generate regular expressions
	#print("{}\t{}".format(seq_filename,seq_key))
	regex_seq1 = replace_N(cur_seq[0:ref_length], trim_end=True)
	regex_seq2 = replace_N(revcomp(cur_seq[-ref_length:]), trim_end=True)
	regex_seqs[seq_key] = [regex_seq1, regex_seq2]

# Output regex references FASTA-like format
regex_file = open(out_filename, "w")
for seq_key in sorted(regex_seqs.keys()):
	regex_file.write(">" + seq_key + "\n")
	regex_file.write(regex_seqs[seq_key][0] + "\n")
	regex_file.write(regex_seqs[seq_key][1] + "\n")
regex_file.close()
