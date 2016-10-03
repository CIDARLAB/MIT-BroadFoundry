#!/usr/bin/env python

#	Copyright (C) 2016 by
#	Alexander Cristofaro <acristo@broadinstitute.org>, Broad Foundry, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

"""
Picard Mark duplicates for each sample
GATK HaplotypeCaller on all samples
"""

# Supporting modules
import argparse
import genetic_analyzer as ga
import sys, os

def main():
    # Parse the command line inputs
    parser = argparse.ArgumentParser(description="mark duplicates and call variants")
    parser.add_argument("-settings", dest="settings", required=True, help="settings.txt", metavar="string")
    parser.add_argument("-samples", dest="samples", required=True, help="1,2", metavar="string")
    args = parser.parse_args()
    # Run the command
    samples = args.samples.split(',')
    settings = ga.load_settings(args.settings)
    for s in samples:
        if not os.path.abspath(settings[s]['fasta_file'].rstrip(".fasta") + '.dict'):
            ga.make_ref_dictionary(settings, s)
        status1 = ga.mark_dupes(settings, s)
    if status1 == 1:
        print("GATK mark duplicates failed")
        return 1

    for s in samples:
        status2 = ga.add_readgroups(settings, s)
    if status2 == 1:
        print("Picard Addreadgroups failed")
        return 1

    for s in samples:
        status3 = ga.gatk_haplotype_caller(settings, s)
    if status3 == 1:
        print("GATK HaplotypeCaller failed")
        return 1

    return status1 + status2 + status3


if __name__ == "__main__":
    status = main()
    sys.exit(status)