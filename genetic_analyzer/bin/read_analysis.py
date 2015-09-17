#!/usr/bin/env python

#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import part_characterizer as pc

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="read_analysis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	args = parser.parse_args()
	# Run the command
	settings = pc.load_settings(args.settings)
	counts = {}
	mapped_reads = {}
	sample_names = settings.keys()
	for s in sample_names:
		counts[s] = pc.read_count_file(pc.count_filename(settings, s))
		mapped_reads[s] = pc.load_mapped_reads(settings, s)
	count_matrix = pc.combine_counts(counts, sample_names)
	pc.save_count_matrix(count_matrix, sample_names, pc.count_matrix_filename(settings))
	pc.save_mapped_reads_matrix(mapped_reads, sample_names, pc.mapped_reads_matrix_filename(settings))

if __name__ == "__main__":
	main()
