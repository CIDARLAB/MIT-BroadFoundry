#!/usr/bin/env python
"""
Clean dialout designs
=====================
	Only required for old data.
"""
from __future__ import print_function, division
import os
import sys
import string
import timeit

__author__  = 'Thomas E. Gorochowski, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

## MAIN =======================================================================

start_time = timeit.default_timer()

# Parse command line parameters
if len(sys.argv) != 4:
	print("Usage: python {} <input prefix> <regex prefix> <output prefix>".format(sys.argv[0]), file=sys.stderr)
	sys.exit()
input_prefix, regex_prefix, output_prefix = sys.argv[1:]
output_prefix = output_prefix.strip()

# -----------------------------------------------------------------------------

barcodes_filename = input_prefix + "_dialout_design_unique_barcodes.csv"
output_filename = output_prefix + "_dialout_design_unique_barcodes.csv"

barcodes = {}
with open(barcodes_filename, "rU") as barcode_file:
	# Ignore header
	barcode_file.next()
	for row in barcode_file:
		row_parts = row.split(",")
		if len(row_parts) >= 3:
			if row_parts[0] not in barcodes.keys():
				barcodes[row_parts[0]] = {}
			if row_parts[1] not in barcodes[row_parts[0]].keys():
				barcodes[row_parts[0]][row_parts[1]] = row_parts[2].strip()

# Output part statistics data
# Barcode, # Designs, Design Names...
print("Rewriting dialout barcode unique designs...", file=sys.stdout)
out_file = open(output_filename, "w")
out_file.write("Design,Barcode,# Reads\n")
for design in barcodes.keys():
	for barcode in barcodes[design].keys():
		out_list = [ design, barcode, barcodes[design][barcode] ]
		out_file.write(",".join(out_list) + "\n")
out_file.close()

# -----------------------------------------------------------------------------

regex_filename = regex_prefix + "_regexs.txt"
designs_filename = input_prefix + "_dialout_designs.csv"
output_filename = output_prefix + "_dialout_designs.csv"

all_designs = []
with open(regex_filename, "rU") as regex_file:
	for design in regex_file:
		design_name = design[1:].strip()
		# Skip regexs
		regex_file.next()
		regex_file.next()
		if design_name != "":
			if design_name not in all_designs:
				all_designs.append(design_name)

designs = {}
with open(designs_filename, "rU") as designs_file:
	# Ignore header
	designs_file.next()
	for row in designs_file:
		row_parts = row.split(",")
		if len(row_parts) >= 3:
			designs[row_parts[0]] = [row_parts[1], row_parts[2].strip()]

# Output part statistics data
# Barcode, # Designs, Design Names...
print("Rewriting dialout designs...", file=sys.stdout)
out_file = open(output_filename, "w")
out_file.write("Design,# Barcodes,# Unique Barcodes\n")
for design in designs.keys():
	out_list = [ design ] + designs[design]
	out_file.write(",".join(out_list) + "\n")
for design in all_designs:
	if design not in designs.keys():
		out_list = [ design ] + ["0", "0"]
		out_file.write(",".join(out_list) + "\n")
out_file.close()

# -----------------------------------------------------------------------------

stop_time = timeit.default_timer()
print("Done ({0:.2f} seconds)".format(stop_time-start_time), file=sys.stdout)
