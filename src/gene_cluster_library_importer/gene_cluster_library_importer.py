#!/usr/bin/env python
"""
Generate a Gene Cluster Library from CSV files
==============================================

    This script will generate a Gene Cluster Library file that contains
    all design and part information. The appropriate spreadsheets must
    be filled in and exported to CSV format for processing by this script.

    Usage:
    ------
    python gene_library_importer.py PART_FILE VARIANT_DESIGNS_FILE 
           VARIANT_ATTRIBS_FILE VARIANT_PARTS_ATTRIBS_FILE OUTPUT_FILE
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import sys
import getopt
import csv
import numpy as np
import gene_cluster_tools as gct

def __make_float_if_needed(s):
	try:
		float(s)
		return float(s)
	except ValueError:
		return s

def load_part_attribs (parts_filename):
	part_data = {}
	parts_reader = csv.reader(open(parts_filename, 'rU'), delimiter=',')
	header = next(parts_reader)
	header_map = {}
	for i in range(len(header)):
		header_map[header[i]] = i
	attrib_keys = [k for k in header_map.keys() if k not in ['Part', 'Part Type', 'Part Sequence']]
	#for k in header_map.keys():
	#	if k not in ['Part', 'Part Type', 'Part Sequence']:
	#		attrib_idxs.append(k)
	for row in parts_reader:
		# Make the attributes map
		part_attribs_map = {}
		for k in attrib_keys:
			part_attribs_map[k] = __make_float_if_needed(row[header_map[k]])
		part_name = row[header_map['Part']]
		part_type = row[header_map['Part Type']]
		part_seq = row[header_map['Part Sequence']]
		part_data[part_name] = [part_type, part_seq, part_attribs_map]
	return part_data

def load_variant_attribs (variant_attribs_filename):
	variant_attribs_data = {}
	variant_attribs_reader = csv.reader(open(variant_attribs_filename, 'rU'), delimiter=',')
	header = next(variant_attribs_reader)
	header_map = {}
	for i in range(len(header)):
		header_map[header[i]] = i
	attrib_keys = [k for k in header_map.keys() if k not in ['Variant']]
	for row in variant_attribs_reader:
		variant_attribs_map = {}
		for k in attrib_keys:
			variant_attribs_map[k] = __make_float_if_needed(row[header_map[k]])
		variant_name = row[header_map['Variant']]
		variant_attribs_data[variant_name] = variant_attribs_map
	return variant_attribs_data

def load_variant_designs (variant_designs_filename):
	# Dictionary keyed on ID and then list returned of all parts present
	variant_designs_data = {}
	variant_designs_reader = csv.reader(open(variant_designs_filename, 'rU'), delimiter=',')
	# Ignore header (variant ID is always in column 0)
	next(variant_designs_reader, None)
	for row in variant_designs_reader:
		variant_designs_data[row[0]] = []
		d = [x.strip() for x in row[1::]]
		for i in d:
			if i != '':
				variant_designs_data[row[0]].append(i)
	return variant_designs_data

def load_variant_part_attribs (variant_part_attribs_filename):
	variant_part_attribs_data = {}
	variant_part_attribs_reader = csv.reader(open(variant_part_attribs_filename, 'rU'), delimiter=',')
	header = next(variant_part_attribs_reader)
	header_map = {}
	for i in range(len(header)):
		header_map[header[i]] = i
	attrib_keys = [k for k in header_map.keys() if k not in ['Variant']]
	for row in variant_part_attribs_reader:
		variant_name = row[header_map['Variant']]
		attrib_name = row[header_map['Attribute']]
		for part in [k for k in header_map.keys() if k not in ['Variant', 'Attribute']]:
			if variant_name not in variant_part_attribs_data.keys():
				variant_part_attribs_data[variant_name] = {}
			if part not in variant_part_attribs_data[variant_name].keys():
				variant_part_attribs_data[variant_name][part] = {}
			variant_part_attribs_data[variant_name][part][attrib_name] = __make_float_if_needed(row[header_map[part]])
	return variant_part_attribs_data

def make_library_from_csv (part_attribs_fn, variant_designs_fn, variant_attribs_fn, variant_part_attribs_fn, output_fn):
	# Load the data
	part_attribs = load_part_attribs(part_attribs_fn)
	variant_designs = load_variant_designs(variant_designs_fn)
	variant_attribs = load_variant_attribs(variant_attribs_fn)
	variant_part_attribs = load_variant_part_attribs(variant_part_attribs_fn)
	# Create GeneClusterLibrary object
	gene_cluster_lib = gct.GeneClusterLibrary()
	# Add the part data
	for p_key in part_attribs.keys():
		gene_cluster_lib.new_part(p_key, part_attribs[p_key][0], part_attribs[p_key][1], attribs=part_attribs[p_key][2])
	# Add the variants
	for v_key in variant_designs.keys():
		part_list = []
		for i in range(len(variant_designs[v_key])):
			cur_part = variant_designs[v_key][i]
			part_el = {}
			if cur_part[0] == 'r':
				part_el['part_name'] = cur_part[1::]
				part_el['dir'] = 'R'
			else:
				part_el['part_name'] = cur_part
				part_el['dir'] = 'F'
			# Load the variant part attributes
			if v_key in variant_part_attribs.keys():
				if part_el['part_name'] in variant_part_attribs[v_key].keys():
					for k in variant_part_attribs[v_key][part_el['part_name']].keys():
						part_el[k] = variant_part_attribs[v_key][part_el['part_name']][k]
			part_list.append(part_el)
		gene_cluster_lib.new_variant(v_key, part_list, seq=None, attribs=variant_attribs[v_key])
	# Save the library
	gene_cluster_lib.save(output_fn)

def main():
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
    # process arguments
    make_library_from_csv(args[0], args[1], args[2], args[3], args[4])

if __name__ == "__main__":
    main()

#make_library_from_csv('part_attributes.csv', 'variant_designs.csv', 'variant_attributes.csv', 'variant_part_attributes.csv', 'test.txt')

