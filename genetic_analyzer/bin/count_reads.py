#!/usr/bin/env python

#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_analyzer as ga
import sys


def main():
    # Parse the command line inputs
    parser = argparse.ArgumentParser(description="count_reads")
    parser.add_argument("-settings", dest="settings", required=True, help="settings.txt", metavar="string")
    parser.add_argument("-samples", dest="samples", required=True, help="1,2", metavar="string")
    parser.add_argument("-feature", dest="feature", required=False, help="gene", metavar="string")
    parser.add_argument("-attribute", dest="attribute", required=True, help="Name", metavar="string")
    parser.add_argument("-strand_opt", dest="strand_opt", required=True, help="no/yes/reverse", metavar="string")
    args = parser.parse_args()
    # Run the command
    samples = args.samples.split(',')
    settings = ga.load_settings(args.settings)
    a = args.attribute
    s_opt = args.strand_opt
    for s in samples:
        features = ga.load_features(settings, s)
        status1 = ga.mapped_reads(settings, s)
        status2 = ga.gene_lengths(settings, s, features=features, attribute=a)
        for f in features:
            status3 = ga.count_reads(settings, s, feature=f, attribute=a, strand_opt=s_opt)
            if status3 == 0:
                continue
            else:
                return 1
    return (status1 + status2)


if __name__ == "__main__":
    status = main()
    sys.exit(status)
