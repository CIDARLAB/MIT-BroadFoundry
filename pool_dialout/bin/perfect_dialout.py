#!/usr/bin/env python
"""
Dialout barcodes from pool
==========================
	Allows for the extraction of barcodes from a pool where design structure is
	fixed. Will only search for perfect matches with references given. Refs are
	given in the form of regular expressions that are matched on the pair of 
	reads. Matched groups should be the barcodes and are recoded against the 
	matched read for reporting. Forward and reverse primer for generating the
	sequencing library are required to orientate the pair of reads.

	Important to consider that barcodes are given as found starting from the 
	forward primer 0...n (barcodes on the paired reads are therefore reversed).
"""
from __future__ import print_function, division
import itertools
import os
import sys
import re
import string
import timeit

__author__  = 'Tarjei Mikkelsen, BTL, Broad Institute\n\
               Thomas E. Gorochowski, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

## HELPERS ====================================================================

# Reverse complement
def revcomp(seq, trans=string.maketrans("ACGT", "TGCA")):
	return "".join(reversed(seq.translate(trans)))

## MAIN =======================================================================

start_time = timeit.default_timer()

# Parse command line parameters
if len(sys.argv) != 9:
	print("Usage: python {} <design regexs> <R1 fastq> <R2 fastq> <fwd primer> <rev primer> <output prefix>".format(sys.argv[0]), file=sys.stderr)
	sys.exit()
design_filename, r1_filename, r2_filename, fwd_primer_len, rev_primer_len, fwd_bc_idx, rev_bc_idx, out_prefix = sys.argv[1:]
out_prefix = out_prefix.strip()

# Allow user to specify index starting at 1 (not 0)
fwd_bc_idx = int(fwd_bc_idx) - 1
rev_bc_idx = int(rev_bc_idx) - 1
fwd_primer_len = int(fwd_primer_len)
rev_primer_len = int(rev_primer_len)

# Load regular expressions uniquely defining each design and compile
design_regexs = {}
fwd_primer = ''
rev_primer = ''
with open(design_filename, "rU") as design_file:
	for header in design_file:
		design_name = header[1:].strip()
		regex_fwd = design_file.next().strip()
		regex_rev = design_file.next().strip()
		design_regexs[design_name] = (re.compile(regex_fwd), re.compile(regex_rev))
		# set the primer
		fwd_primer = regex_fwd[0:fwd_primer_len]
		rev_primer = regex_rev[0:rev_primer_len]

# Process the raw reads and try to match
print("Starting to process read...", file=sys.stdout)
sys.stdout.flush()
n_reads = 0
n_accepted = 0
n_matched = 0
found_designs = {}
found_barcodes = {}
found_barcode_reads = {}
barcodes_to_check = [{},{}]
with open(r1_filename, "rU") as r1, open(r2_filename, "rU") as r2:
	r1_r2 = itertools.izip(r1, r2)
	for header1, header2 in r1_r2:

		# Give indication of progress
		if n_reads % 1000 == 0:
			print("Processed {} reads".format(n_reads), file=sys.stdout)
			sys.stdout.flush()

		# Extract data and clean
		seq1, seq2 = r1_r2.next()
		plus1, plus2 = r1_r2.next()
		qual1, qual2 = r1_r2.next()
		seq1, seq2 = seq1.strip(), seq2.strip()
		qual1, qual2 = qual1.strip(), qual2.strip()
		# Check that paired-end read
		read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
		assert read_name1 == read_name2
		n_reads += 1
		# Flip reads if the reverse primer is seen read first
		if (seq1.startswith(rev_primer) and
			seq2.startswith(fwd_primer)):
			seq1, seq2 = seq2, seq1
			qual1, qual2 = qual2, qual1
		# Reject reads if we don't see perfect primer matches
		if not (seq1.startswith(fwd_primer) and
				seq2.startswith(rev_primer)):
			continue
		n_accepted += 1

		# Attempt to match to regexs
		found_design = None
		found_barcode = None
		for design in design_regexs.keys():
			m1 = re.match(design_regexs[design][0], seq1)
			m2 = re.match(design_regexs[design][1], seq2)
			# Assumes paired matches are unique
			if (m1 != None) and (m2 != None):
				found_design = design
				found_barcode = list(m1.groups()) + [revcomp(x) for x in m2.groups()]
				barcode_key = "-".join(found_barcode)
				
				# CHECK BARCODES ##################################################
				if found_barcode[fwd_bc_idx] not in barcodes_to_check[0].keys():
					barcodes_to_check[0][found_barcode[fwd_bc_idx]] = [found_design]
				else:
					if found_design not in barcodes_to_check[0][found_barcode[fwd_bc_idx]]:
						barcodes_to_check[0][found_barcode[fwd_bc_idx]].append(found_design)
				if found_barcode[rev_bc_idx] not in barcodes_to_check[1].keys():
					barcodes_to_check[1][found_barcode[rev_bc_idx]] = [found_design]
				else:
					if found_design not in barcodes_to_check[1][found_barcode[rev_bc_idx]]:
						barcodes_to_check[1][found_barcode[rev_bc_idx]].append(found_design)
				###############################################################

				if barcode_key not in found_barcodes.keys():
					found_barcodes[barcode_key] = {}
					found_barcodes[barcode_key][found_design] = 1
				else:
					if found_design in found_barcodes[barcode_key].keys():
						found_barcodes[barcode_key][found_design] += 1
					else:
						found_barcodes[barcode_key][found_design] = 1
				break
		# If found match then process else continue to next read
		if found_design is None:
			continue
		n_matched += 1

		# Only add barcode if not already seen
		if found_design not in found_designs.keys():
			found_designs[found_design] = [found_barcode]
		else:
			if found_barcode not in found_designs[found_design]:
				found_designs[found_design].append(found_barcode)

# Output the matched designs and barcode data
print("Writing dialout for designs...", file=sys.stdout)
# Design, # Reads, # Barcodes, # Unique Barcodes, Unique Barcodes... 
design_file = open(out_prefix + "dialout_designs.csv", "w")
design_file.write("Design,# Barcode Reads,# Unique Barcode Reads\n")
# Design, Barcode, # Reads 
unique_bc_file = open(out_prefix + "dialout_design_unique_barcodes.csv", "w")
unique_bc_file.write("Design,Unique Barcode,# Reads\n")
n_unique_bc = 0
design_with_unique_bc = []
for design in sorted(found_designs.keys()):
	print("Processing design: {}".format(design), file=sys.stdout)
	sys.stdout.flush()
	unique_barcodes = []
	unique_for_design = False
	for bc in found_designs[design]:
		# Check to see where barcode is used (if unique for design)
		keys_0 = barcodes_to_check[0].keys()
		keys_1 = barcodes_to_check[1].keys()
		if ( (len(barcodes_to_check[0][bc[fwd_bc_idx]]) == 1) and
		     (revcomp(bc[fwd_bc_idx]) not in keys_0) and
		     (bc[fwd_bc_idx] not in keys_1) and
		     (revcomp(bc[fwd_bc_idx]) not in keys_1) and
		     (len(barcodes_to_check[1][bc[rev_bc_idx]]) == 1) and
		     (revcomp(bc[rev_bc_idx]) not in keys_1) and
		     (bc[rev_bc_idx] not in keys_0) and
		     (revcomp(bc[rev_bc_idx]) not in keys_0) ):
			unique_barcodes.append(bc)
			n_unique_bc += 1
			unique_for_design = True
	# For summary statistic of designs with > 0 unique barcodes
	if unique_for_design == True:
		design_with_unique_bc.append(design)
	out_list = [design, str(len(found_designs[design])), str(len(unique_barcodes))]
	# For the unique barcode add to the list or output
	for bc in unique_barcodes:
		unique_bc_file.write(",".join([design, "-".join(bc), str(found_barcodes["-".join(bc)][design])]) + "\n") 
	design_file.write(",".join(out_list) + "\n")
# Add designs not matched at end
for design in design_regexs.keys():
	if design not in found_designs.keys():
		out_list = [ design ] + ["0", "0"]
		design_file.write(",".join(out_list) + "\n")
design_file.close()
unique_bc_file.close()

# Output barcode data
# Barcode, # Designs, Design Names...
print("Writing dialout for barcodes...", file=sys.stdout)
barcode_file = open(out_prefix + "dialout_barcodes.csv", "w")
barcode_file.write("Barcode ID,# Designs,Design Names\n")
for barcode_id in sorted(found_barcodes.keys()):
	out_list = [barcode_id, str(len(found_barcodes[barcode_id]))] + found_barcodes[barcode_id].keys()
	barcode_file.write(",".join(out_list) + "\n")
barcode_file.close()

# Output summary statistics to log
print("Writing summary...", file=sys.stdout)
log = open(out_prefix + "dialout_summary.txt", "w")
print("READ SUMMARY", file=log)
print("Total read pairs: \t{}".format(n_reads), file=log)
print("Valid read pairs: \t{}".format(n_accepted), file=log)
print("Matched read pairs: \t{}".format(n_matched), file=log)
print("", file=log)
print("COVERAGE SUMMARY", file=log)
print("Matched reads: \t{0:.2f}%".format((float(n_matched)/n_reads)*100.0), file=log)
print("Unique barcodes: \t{0:.2f}%".format((float(n_unique_bc)/n_matched)*100.0), file=log)
print("Uniquely barcoded designs: \t{0:.2f}%".format((float(len(design_with_unique_bc))/len(design_regexs.keys()))*100.0), file=log)
print("", file=log)
print("RUNTIME", file=log)
stop_time = timeit.default_timer()
print("Total time: \t{0:.2f}%".format(stop_time-start_time), file=log)
log.close()

print("Done ({0:.2f} seconds)".format(stop_time-start_time), file=sys.stdout)
print("")
sys.stdout.flush()

