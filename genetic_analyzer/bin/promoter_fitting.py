#!/usr/bin/env python

#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_analyzer as ga

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="map_reads")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-output_name",  dest="output_name",  required=True,  help="sample_set1", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	args = parser.parse_args()
	# Run the command
	samples = args.samples.split(',')
	settings = ga.load_settings(args.settings)
	ga.fit_promoter_response_functions(settings, samples, args.output_name)

if __name__ == "__main__":
	main()
