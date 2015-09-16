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
	parser = argparse.ArgumentParser(description="count_reads")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	parser.add_argument("-feature",  dest="feature",  required=True,  help="gene", metavar="string")
	parser.add_argument("-attribute",  dest="attribute",  required=True,  help="Name", metavar="string")
	args = parser.parse_args()
	# Run the command
	samples = args.samples.split(',')
	settings = pc.load_settings(args.settings)
	f = args.feature
	a = args.attribute
	for s in samples:
		pc.count_reads(settings, s, feature=f, attribute=a)
		pc.gene_lengths(settings, s, feature=f, attribute=a)

if __name__ == "__main__":
	main()
