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
import gene_cluster_tools as gct

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
    gct.make_library_from_csv(args[0], args[1], args[2], args[3], args[4])

if __name__ == "__main__":
    main()

#make_library_from_csv('part_attributes.csv', 'variant_designs.csv', 'variant_attributes.csv', 'variant_part_attributes.csv', 'test.txt')

