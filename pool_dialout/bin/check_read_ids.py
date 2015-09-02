#!/usr/bin/env python
"""
Check the names of pairs of FASTQ files
"""
from __future__ import print_function, division
import itertools
import os
import sys
import re
import string
import timeit

__author__  = 'Thomas E. Gorochowski, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

## MAIN =======================================================================

start_time = timeit.default_timer()

# Parse command line parameters
if len(sys.argv) != 3:
	print("Usage: python {} <R1 fastq> <R2 fastq>".format(sys.argv[0]), file=sys.stderr)
	sys.exit()
r1_filename, r2_filename = sys.argv[1:]
r1_filename = r1_filename.strip()
r2_filename = r2_filename.strip()

# Read it all in (hope there is enough memory)
file_r1 = open(r1_filename, "rU")
file_r2 = open(r2_filename, "rU")
print("Reading file 1...")
r1_content = file_r1.readlines()
print("Reading file 2...")
r2_content = file_r2.readlines()

print("")
print("R1 length (lines) = {}".format(len(r1_content)))
print("R2 length (lines) = {}".format(len(r2_content)))
print("")
sys.stdout.flush()

line_idx = 0
max_line_idx = len(r1_content)

while line_idx < max_line_idx:
	# Extract data and clean
	header1 = r1_content[line_idx]
	header2 = r2_content[line_idx]
	seq1 = r1_content[line_idx+1]
	seq2 = r2_content[line_idx+1]
	plus1 = r1_content[line_idx+2]
	plus2 = r2_content[line_idx+2]
	qual1 = r1_content[line_idx+3]
	qual2 = r2_content[line_idx+3]
	line_idx += 4
	seq1, seq2 = seq1.strip(), seq2.strip()
	qual1, qual2 = qual1.strip(), qual2.strip()
	# Check that paired-end read
	read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
	
	if read_name1 != read_name2:
		print("ERROR at line {}".format(line_idx))
		print("   R1: {}".format(read_name1))
		print("       {}".format(header1))
		print("       {}".format(seq1))
		print("       {}".format(plus1))
		print("       {}".format(qual1))
		print("   R2: {}".format(read_name2))
		print("       {}".format(header2))
		print("       {}".format(seq2))
		print("       {}".format(plus2))
		print("       {}".format(qual2))
		print("")
		sys.stdout.flush()

sys.stdout.flush()
