#!/usr/bin/env python
"""
Generate a GeneClusterLibrary for the circuit and clean the trace
files to make them compatible with other scripts. We hard code everything
as it is easier and changes can still be incorporated fairly easily, e.g.,
performance data for the parts.

IMPORTANT: we reformat the sequencing results so that the last 1000bp are
concatentated to the start of the construct. This is to enable easier
analysis of promoters on the edge of the sequencing data. This doesn't
cause any issues as this is a circular plasmid after all.

For reference the naming of samples is:
	t1 = none
	t2 = iptg                 (PTac)
	t3 = atc                  (PTet)
	t4 = iptg + atc           (PTac, PTet)
	t5 = ara                  (PBAD)
	t6 = ara + iptg           (PBAD, PTac)
	t7 = ara + atc            (PBAD, PTet)
	t8 = ara + iptg + atc     (PBAD, PTac, PTet)

In addition, REU values for the working rep 1 and non-working rep 2 are
included as attributes. Specific parts for the repressors are:
	AmtR rbs1
	LitR rbs2
	BM3R1 rbs1
	SrpR rbs0
	PhlF rbs1
	YFP rbs0
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import gene_cluster_library as gcl

import numpy as np
import csv

# Create GeneClusterLibrary object
circuit_lib = gcl.GeneClusterLibrary()

######################################################################
# ADD THE PARTS (we only consider transcriptional elements)
######################################################################

# Promoters
circuit_lib.new_part('PTac',   'Promoter', 'TGTTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGCTCACAATT', attribs=None)
circuit_lib.new_part('PTet',   'Promoter', 'TTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCAC', attribs=None)
circuit_lib.new_part('PBAD',   'Promoter', 'AGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGGCTAGC', attribs=None)
circuit_lib.new_part('PBM3R1', 'Promoter', 'TCTGATTCGTTACCAATTGACGGAATGAACGTTCATTCCGATAATGCTAGC', attribs=None)
circuit_lib.new_part('PAmtR',  'Promoter', 'GATTCGTTACCAATTGACAGTTTCTATCGATCTATAGATAATGCTAGC', attribs=None)
circuit_lib.new_part('PSrpR',  'Promoter', 'GATTCGTTACCAATTGACAGCTAGCTCAGTCCTAGGTATATACATACATGCTTGTTTGTTTGTAAAC', attribs=None)
circuit_lib.new_part('PLitR',  'Promoter', 'GATTCGTTACCAATTGACAAATTTATAAATTGTCAGTATAATGCTAGC', attribs=None)
circuit_lib.new_part('PPhlF',  'Promoter', 'TCTGATTCGTTACCAATTGACATGATACGAAACGTACCGTATCGTTAAGGT', attribs=None)

# CDSs
circuit_lib.new_part('AmtR',  'CDS', 'ATGGCAGGCGCAGTTGGTCGTCCGCGTCGTAGTGCACCGCGTCGTGCAGGTAAAAATCCGCGTGAAGAAATTCTGGATGCAAGCGCAGAACTGTTTACCCGTCAGGGTTTTGCAACCACCAGTACCCATCAGATTGCAGATGCAGTTGGTATTCGTCAGGCAAGCCTGTATTATCATTTTCCGAGCAAAACCGAAATCTTTCTGACCCTGCTGAAAAGCACCGTTGAACCGAGCACCGTTCTGGCAGAAGATCTGAGCACCCTGGATGCAGGTCCGGAAATGCGTCTGTGGGCAATTGTTGCAAGCGAAGTTCGTCTGCTGCTGAGCACCAAATGGAATGTTGGTCGTCTGTATCAGCTGCCGATTGTTGGTAGCGAAGAATTTGCAGAATATCATAGCCAGCGTGAAGCACTGACCAATGTTTTTCGTGATCTGGCAACCGAAATTGTTGGTGATGATCCGCGTGCAGAACTGCCGTTTCATATTACCATGAGCGTTATTGAAATGCGTCGCAATGATGGTAAAATTCCGAGTCCGCTGAGCGCAGATAGCCTGCCGGAAACCGCAATTATGCTGGCAGATGCAAGCCTGGCAGTTCTGGGTGCACCGCTGCCTGCAGATCGTGTTGAAAAAACCCTGGAACTGATTAAACAGGCAGATGCAAAATAA', attribs=None)
circuit_lib.new_part('LitR',  'CDS', 'ATGGATACCATTCAGAAACGTCCGCGTACCCGTCTGAGTCCGGAAAAACGTAAAGAACAGCTGCTGGATATTGCCATTGAAGTTTTTAGCCAGCGTGGTATTGGTCGTGGTGGTCATGCAGATATTGCAGAAATTGCACAGGTTAGCGTTGCAACCGTGTTTAACTATTTTCCGACCCGTGAAGATCTGGTTGATGATGTTCTGAACAAAGTGGAAAACGAGTTTCACCAGTTCATCAATAACAGCATTAGCCTGGATCTGGATGTTCGTAGCAATCTGAATACCCTGCTGCTGAACATTATTGATAGCGTTCAGACCGGCAACAAATGGATTAAAGTTTGGTTTGAATGGTCAACCAGCACCCGTGATGAAGTTTGGCCTCTGTTTCTGAGCACCCATAGCAATACCAATCAGGTGATCAAAACCATGTTTGAAGAGGGTATTGAACGCAATGAAGTGTGCAATGATCATACACCGGAAAATCTGACCAAAATGCTGCATGGTATTTGCTATAGCGTGTTTATTCAGGCCAATCGTAATAGCAGCAGCGAAGAAATGGAAGAAACCGCAAATTGCTTTCTGAATATGCTGTGCATCTACAAATAA', attribs=None)
circuit_lib.new_part('BM3R1', 'CDS', 'ATGGAAAGCACCCCGACCAAACAGAAAGCAATTTTTAGCGCAAGCCTGCTGCTGTTTGCAGAACGTGGTTTTGATGCAACCACCATGCCGATGATTGCAGAAAATGCAAAAGTTGGTGCAGGCACCATTTATCGCTATTTCAAAAACAAAGAAAGCCTGGTGAACGAACTGTTTCAGCAGCATGTTAATGAATTTCTGCAGTGTATTGAAAGCGGTCTGGCAAATGAACGTGATGGTTATCGTGATGGCTTTCATCACATTTTTGAAGGTATGGTGACCTTTACCAAAAATCATCCGCGTGCACTGGGTTTTATCAAAACCCATAGCCAGGGCACCTTTCTGACCGAAGAAAGCCGTCTGGCATATCAGAAACTGGTTGAATTTGTGTGCACCTTTTTTCGTGAAGGTCAGAAACAGGGTGTGATTCGTAATCTGCCGGAAAATGCACTGATTGCAATTCTGTTTGGCAGCTTTATGGAAGTGTATGAAATGATCGAGAACGATTATCTGAGCCTGACCGATGAACTGCTGACCGGTGTTGAAGAAAGCCTGTGGGCAGCACTGAGCCGTCAGAGCTAA', attribs=None)
circuit_lib.new_part('SrpR',  'CDS', 'ATGGCACGTAAAACCGCAGCAGAAGCAGAAGAAACCCGTCAGCGTATTATTGATGCAGCACTGGAAGTTTTTGTTGCACAGGGTGTTAGTGATGCAACCCTGGATCAGATTGCACGTAAAGCCGGTGTTACCCGTGGTGCAGTTTATTGGCATTTTAATGGTAAACTGGAAGTTCTGCAGGCAGTTCTGGCAAGCCGTCAGCATCCGCTGGAACTGGATTTTACACCGGATCTGGGTATTGAACGTAGCTGGGAAGCAGTTGTTGTTGCAATGCTGGATGCAGTTCATAGTCCGCAGAGCAAACAGTTTAGCGAAATTCTGATTTATCAGGGTCTGGATGAAAGCGGTCTGATTCATAATCGTATGGTTCAGGCAAGCGATCGTTTTCTGCAGTATATTCATCAGGTTCTGCGTCATGCAGTTACCCAGGGTGAACTGCCGATTAATCTGGATCTGCAGACCAGCATTGGTGTTTTTAAAGGTCTGATTACCGGTCTGCTGTATGAAGGTCTGCGTAGCAAAGATCAGCAGGCACAGATTATCAAAGTTGCACTGGGTAGCTTTTGGGCACTGCTGCGTGAACCGCCTCGTTTTCTGCTGTGTGAAGAAGCACAGATTAAACAGGTGAAATCCTTCGAATAA', attribs=None)
circuit_lib.new_part('PhlF',  'CDS', 'ATGGCACGTACCCCGAGCCGTAGCAGCATTGGTAGCCTGCGTAGTCCGCATACCCATAAAGCAATTCTGACCAGCACCATTGAAATCCTGAAAGAATGTGGTTATAGCGGTCTGAGCATTGAAAGCGTTGCACGTCGTGCCGGTGCAAGCAAACCGACCATTTATCGTTGGTGGACCAATAAAGCAGCACTGATTGCCGAAGTGTATGAAAATGAAAGCGAACAGGTGCGTAAATTTCCGGATCTGGGTAGCTTTAAAGCCGATCTGGATTTTCTGCTGCGTAATCTGTGGAAAGTTTGGCGTGAAACCATTTGTGGTGAAGCATTTCGTTGTGTTATTGCAGAAGCACAGCTGGACCCTGCAACCCTGACCCAGCTGAAAGATCAGTTTATGGAACGTCGTCGTGAGATGCCGAAAAAACTGGTTGAAAATGCCATTAGCAATGGTGAACTGCCGAAAGATACCAATCGTGAACTGCTGCTGGATATGATTTTTGGTTTTTGTTGGTATCGCCTGCTGACCGAACAGCTGACCGTTGAACAGGATATTGAAGAATTTACCTTCCTGCTGATTAATGGTGTTTGTCCGGGTACACAGCGTTAA', attribs=None)
circuit_lib.new_part('YFP',   'CDS', 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA', attribs=None)

# Terminators
circuit_lib.new_part('L3S2P55',      'Terminator', 'CTCGGTACCAAAGACGAACAATAAGACGCTGAAAAGCGTCTTTTTTCGTTTTGGTCC', attribs=None)
circuit_lib.new_part('L3S2P24',      'Terminator', 'CTCGGTACCAAATTCCAGAAAAGACACCCGAAAGGGTGTTTTTTCGTTTTGGTCC', attribs=None)
circuit_lib.new_part('L3S2P11',      'Terminator', 'CTCGGTACCAAATTCCAGAAAAGAGACGCTTTCGAGCGTCTTTTTTCGTTTTGGTCC', attribs=None)
circuit_lib.new_part('ECK120029600', 'Terminator', 'TTCAGCCAAAAAACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTCTTCTCAA', attribs=None)
circuit_lib.new_part('ECK120033737', 'Terminator', 'GGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGG', attribs=None)
circuit_lib.new_part('L3S2P21',      'Terminator', 'CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC', attribs=None)

# Ribozymes (insulators)
circuit_lib.new_part('BydvJ', 'Ribozyme', 'AGGGTGTCTCAAGGTGCGTACCTTGACTGATGAGTCCGAAAGGACGAAACACCCCTCTACAAATAATTTTGTTTAA', attribs=None)
circuit_lib.new_part('PlmJ', 'Ribozyme', 'AGTCATAAGTCTGGGCTAAGCCCACTGATGAGTCGCTGAAATGCGACGAAACTTATGACCTCTACAAATAATTTTGTTTAA', attribs=None)
circuit_lib.new_part('SarJ', 'Ribozyme', 'AGACTGTCGCCGGATGTGTATCCGACCTGACGATGGCCCAAAAGGGCCGAAACAGTCCTCTACAAATAATTTTGTTTAA', attribs=None)
circuit_lib.new_part('RiboJ10', 'Ribozyme', 'AGCGCTCAACGGGTGTGCTTCCCGTTCTGATGAGTCCGTGAGGACGAAAGCGCCTCTACAAATAATTTTGTTTAA', attribs=None)
circuit_lib.new_part('RiboJ53', 'Ribozyme', 'AGCGGTCAACGCATGTGCTTTGCGTTCTGATGAGACAGTGATGTCGAAACCGCCTCTACAAATAATTTTGTTTAA', attribs=None)
circuit_lib.new_part('RiboJ', 'Ribozyme', 'AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA', attribs=None)



######################################################################
# MAKE THE CIRCUIT (Include FPKM data)
######################################################################

prefix_len = 1000

# Vector sequence
seq = 'GCTTAACGATCGTTGGCTGTGTTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGCTCACAATTTACTCCACCGTTGGCTTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCACCTGAAGGGTGTCTCAAGGTGCGTACCTTGACTGATGAGTCCGAAAGGACGAAACACCCCTCTACAAATAATTTTGTTTAAAATGTTCCCTAATAATCAGCAAAGAGGTTACTAGATGGCAGGCGCAGTTGGTCGTCCGCGTCGTAGTGCACCGCGTCGTGCAGGTAAAAATCCGCGTGAAGAAATTCTGGATGCAAGCGCAGAACTGTTTACCCGTCAGGGTTTTGCAACCACCAGTACCCATCAGATTGCAGATGCAGTTGGTATTCGTCAGGCAAGCCTGTATTATCATTTTCCGAGCAAAACCGAAATCTTTCTGACCCTGCTGAAAAGCACCGTTGAACCGAGCACCGTTCTGGCAGAAGATCTGAGCACCCTGGATGCAGGTCCGGAAATGCGTCTGTGGGCAATTGTTGCAAGCGAAGTTCGTCTGCTGCTGAGCACCAAATGGAATGTTGGTCGTCTGTATCAGCTGCCGATTGTTGGTAGCGAAGAATTTGCAGAATATCATAGCCAGCGTGAAGCACTGACCAATGTTTTTCGTGATCTGGCAACCGAAATTGTTGGTGATGATCCGCGTGCAGAACTGCCGTTTCATATTACCATGAGCGTTATTGAAATGCGTCGCAATGATGGTAAAATTCCGAGTCCGCTGAGCGCAGATAGCCTGCCGGAAACCGCAATTATGCTGGCAGATGCAAGCCTGGCAGTTCTGGGTGCACCGCTGCCTGCAGATCGTGTTGAAAAAACCCTGGAACTGATTAAACAGGCAGATGCAAAATAACTCGGTACCAAAGACGAACAATAAGACGCTGAAAAGCGTCTTTTTTCGTTTTGGTCCAATGACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGGCTAGCTACTCCACCGTTGGCTTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCACCTGAAGTCATAAGTCTGGGCTAAGCCCACTGATGAGTCGCTGAAATGCGACGAAACTTATGACCTCTACAAATAATTTTGTTTAAGTCCTATGGACTTTTTCATACAGGAGAACCCTCGATGGATACCATTCAGAAACGTCCGCGTACCCGTCTGAGTCCGGAAAAACGTAAAGAACAGCTGCTGGATATTGCCATTGAAGTTTTTAGCCAGCGTGGTATTGGTCGTGGTGGTCATGCAGATATTGCAGAAATTGCACAGGTTAGCGTTGCAACCGTGTTTAACTATTTTCCGACCCGTGAAGATCTGGTTGATGATGTTCTGAACAAAGTGGAAAACGAGTTTCACCAGTTCATCAATAACAGCATTAGCCTGGATCTGGATGTTCGTAGCAATCTGAATACCCTGCTGCTGAACATTATTGATAGCGTTCAGACCGGCAACAAATGGATTAAAGTTTGGTTTGAATGGTCAACCAGCACCCGTGATGAAGTTTGGCCTCTGTTTCTGAGCACCCATAGCAATACCAATCAGGTGATCAAAACCATGTTTGAAGAGGGTATTGAACGCAATGAAGTGTGCAATGATCATACACCGGAAAATCTGACCAAAATGCTGCATGGTATTTGCTATAGCGTGTTTATTCAGGCCAATCGTAATAGCAGCAGCGAAGAAATGGAAGAAACCGCAAATTGCTTTCTGAATATGCTGTGCATCTACAAATAACTCGGTACCAAATTCCAGAAAAGACACCCGAAAGGGTGTTTTTTCGTTTTGGTCCTGTCACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGGCTAGCCTGAAGACTGTCGCCGGATGTGTATCCGACCTGACGATGGCCCAAAAGGGCCGAAACAGTCCTCTACAAATAATTTTGTTTAACTATGGACTATGTTTTTCAAAGACGAAAAACTACTAGATGGAAAGCACCCCGACCAAACAGAAAGCAATTTTTAGCGCAAGCCTGCTGCTGTTTGCAGAACGTGGTTTTGATGCAACCACCATGCCGATGATTGCAGAAAATGCAAAAGTTGGTGCAGGCACCATTTATCGCTATTTCAAAAACAAAGAAAGCCTGGTGAACGAACTGTTTCAGCAGCATGTTAATGAATTTCTGCAGTGTATTGAAAGCGGTCTGGCAAATGAACGTGATGGTTATCGTGATGGCTTTCATCACATTTTTGAAGGTATGGTGACCTTTACCAAAAATCATCCGCGTGCACTGGGTTTTATCAAAACCCATAGCCAGGGCACCTTTCTGACCGAAGAAAGCCGTCTGGCATATCAGAAACTGGTTGAATTTGTGTGCACCTTTTTTCGTGAAGGTCAGAAACAGGGTGTGATTCGTAATCTGCCGGAAAATGCACTGATTGCAATTCTGTTTGGCAGCTTTATGGAAGTGTATGAAATGATCGAGAACGATTATCTGAGCCTGACCGATGAACTGCTGACCGGTGTTGAAGAAAGCCTGTGGGCAGCACTGAGCCGTCAGAGCTAACTCGGTACCAAATTCCAGAAAAGAGACGCTTTCGAGCGTCTTTTTTCGTTTTGGTCCTCTGAATCCGCGTGATAGGTCTGATTCGTTACCAATTGACGGAATGAACGTTCATTCCGATAATGCTAGCCTTGTCCAACCAAATGATTCGTTACCAATTGACAGTTTCTATCGATCTATAGATAATGCTAGCCTGAAGCGCTCAACGGGTGTGCTTCCCGTTCTGATGAGTCCGTGAGGACGAAAGCGCCTCTACAAATAATTTTGTTTAACTATGGACTATGTTTTCACACAGGAAATACCAGGATGGCACGTAAAACCGCAGCAGAAGCAGAAGAAACCCGTCAGCGTATTATTGATGCAGCACTGGAAGTTTTTGTTGCACAGGGTGTTAGTGATGCAACCCTGGATCAGATTGCACGTAAAGCCGGTGTTACCCGTGGTGCAGTTTATTGGCATTTTAATGGTAAACTGGAAGTTCTGCAGGCAGTTCTGGCAAGCCGTCAGCATCCGCTGGAACTGGATTTTACACCGGATCTGGGTATTGAACGTAGCTGGGAAGCAGTTGTTGTTGCAATGCTGGATGCAGTTCATAGTCCGCAGAGCAAACAGTTTAGCGAAATTCTGATTTATCAGGGTCTGGATGAAAGCGGTCTGATTCATAATCGTATGGTTCAGGCAAGCGATCGTTTTCTGCAGTATATTCATCAGGTTCTGCGTCATGCAGTTACCCAGGGTGAACTGCCGATTAATCTGGATCTGCAGACCAGCATTGGTGTTTTTAAAGGTCTGATTACCGGTCTGCTGTATGAAGGTCTGCGTAGCAAAGATCAGCAGGCACAGATTATCAAAGTTGCACTGGGTAGCTTTTGGGCACTGCTGCGTGAACCGCCTCGTTTTCTGCTGTGTGAAGAAGCACAGATTAAACAGGTGAAATCCTTCGAATAATTCAGCCAAAAAACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTCTTCTCAAGGAGTCTATGATTGGTCCAGATTCGTTACCAATTGACAGCTAGCTCAGTCCTAGGTATATACATACATGCTTGTTTGTTTGTAAACCGAGCGTAGAGCTTAGATTCGTTACCAATTGACAAATTTATAAATTGTCAGTATAATGCTAGCCTGAAGCGGTCAACGCATGTGCTTTGCGTTCTGATGAGACAGTGATGTCGAAACCGCCTCTACAAATAATTTTGTTTAAGGAGCTATGGACTATGTTTGAAAGGCTGAAATACTAGATGGCACGTACCCCGAGCCGTAGCAGCATTGGTAGCCTGCGTAGTCCGCATACCCATAAAGCAATTCTGACCAGCACCATTGAAATCCTGAAAGAATGTGGTTATAGCGGTCTGAGCATTGAAAGCGTTGCACGTCGTGCCGGTGCAAGCAAACCGACCATTTATCGTTGGTGGACCAATAAAGCAGCACTGATTGCCGAAGTGTATGAAAATGAAAGCGAACAGGTGCGTAAATTTCCGGATCTGGGTAGCTTTAAAGCCGATCTGGATTTTCTGCTGCGTAATCTGTGGAAAGTTTGGCGTGAAACCATTTGTGGTGAAGCATTTCGTTGTGTTATTGCAGAAGCACAGCTGGACCCTGCAACCCTGACCCAGCTGAAAGATCAGTTTATGGAACGTCGTCGTGAGATGCCGAAAAAACTGGTTGAAAATGCCATTAGCAATGGTGAACTGCCGAAAGATACCAATCGTGAACTGCTGCTGGATATGATTTTTGGTTTTTGTTGGTATCGCCTGCTGACCGAACAGCTGACCGTTGAACAGGATATTGAAGAATTTACCTTCCTGCTGATTAATGGTGTTTGTCCGGGTACACAGCGTTAAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTACGCGACGTACGGTGGAATCTGATTCGTTACCAATTGACATGATACGAAACGTACCGTATCGTTAAGGTAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAATACTAGAGAAAGAGGGGAAATACTAGATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCAATGGTCACCATATATCAAGTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCTACTAGAGAAAGAGGAGAAATACTAGATGGCTGAAGCGCAAAATGATCCCCTGCTGCCGGGATACTCGTTTAATGCCCATCTGGTGGCGGGTTTAACGCCGATTGAGGCCAACGGTTATCTCGATTTTTTTATCGACCGACCGCTGGGAATGAAAGGTTATATTCTCAATCTCACCATTCGCGGTCAGGGGGTGGTGAAAAATCAGGGACGAGAATTTGTTTGCCGACCGGGTGATATTTTGCTGTTCCCGCCAGGAGAGATTCATCACTACGGTCGTCATCCGGAGGCTCGCGAATGGTATCACCAGTGGGTTTACTTTCGTCCGCGCGCCTACTGGCATGAATGGCTTAACTGGCCGTCAATATTTGCCAATACGGGGTTCTTTCGCCCGGATGAAGCGCACCAGCCGCATTTCAGCGACCTGTTTGGGCAAATCATTAACGCCGGGCAAGGGGAAGGGCGCTATTCGGAGCTGCTGGCGATAAATCTGCTTGAGCAATTGTTACTGCGGCGCATGGAAGCGATTAACGAGTCGCTCCATCCACCGATGGATAATCGGGTACGCGAGGCTTGTCAGTACATCAGCGATCACCTGGCAGACAGCAATTTTGATATCGCCAGCGTCGCACAGCATGTTTGCTTGTCGCCGTCGCGTCTGTCACATCTTTTCCGCCAGCAGTTAGGGATTAGCGTCTTAAGCTGGCGCGAGGACCAACGTATCAGCCAGGCGAAGCTGCTTTTGAGCACCACCCGGATGCCTATCGCCACCGTCGGTCGCAATGTTGGTTTTGACGATCAACTCTATTTCTCGCGGGTATTTAAAAAATGCACCGGGGCCAGCCCGAGCGAGTTCCGTGCCGGTTGTGAAGAAAAAGTGAATGATGTAGCCGTCAAGTTGTCATAATAACCAATTATTGAAGGCCGCTAACGCGGCCTTTTTTTGTTTCTGGTCTCCCAATGGCGGCGCGCCATCGAATGGCGCAAAACCTTTCGCGGTATGGCATGATAGCGCCCGGAAGAGAGTCAATTCAGGGTGGTGAATATGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCGGCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAGTCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTCGCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAACGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGCACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATTTTCTCCCATGAGGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAGCAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGCTGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGGAGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACTGCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCCGGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGATAGCTCATGTTATATCCCGCCGTTAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCAGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGATAATCCAGGAGGAAAAAAATGTCCAGATTAGATAAAAGTAAAGTGATTAACAGCGCATTAGAGCTGCTTAATGAGGTCGGAATCGAAGGTTTAACAACCCGTAAACTCGCCCAGAAGCTAGGTGTAGAGCAGCCTACATTGTATTGGCATGTAAAAAATAAGCGGGCTTTGCTCGACGCCTTAGCCATTGAGATGTTAGATAGGCACCATACTCACTTTTGCCCTTTAGAAGGGGAAAGCTGGCAAGATTTTTTACGTAATAACGCTAAAAGTTTTAGATGTGCTTTACTAAGTCATCGCGATGGAGCAAAAGTACATTTAGGTACACGGCCTACAGAAAAACAGTATGAAACTCTCGAAAATCAATTAGCCTTTTTATGCCAACAAGGTTTTTCACTAGAGAATGCATTATATGCACTCAGCGCTGTGGGGCATTTTACTTTAGGTTGCGTATTGGAAGATCAAGAGCATCAAGTCGCTAAAGAAGAAAGGGAAACACCTACTACTGATAGTATGCCGCCATTATTACGACAAGCTATCGAATTATTTGATCACCAAGGTGCAGAGCCAGCCTTCTTATTCGGCCTTGAATTGATCATATGCGGATTAGAAAAACAACTTAAATGTGAAAGTGGGTCCTAATAATTGGTAACGAATCAGACAATTGACGGCTCGAGGGAGTAGCATAGGGTTTGCAGAATCCCTGCTTCGTCCATTTGACAGGCACATTATGCATCGATGATAAGCTGTCAAACATGAGCAGATCCTCTACGCCGGACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATCGCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTCATGAGCAAATATTTTATCTGAGGTGCTTCCTCGCTCACTGACTCGCTGCACGAGGCAGACCTCAGCGCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGAAGTGCTTCATGTGGCAGGAGAAAAAAGGCTGCACCGGTGCGTCAGCAGAATATGTGATACAGGATATATTCCGCTTCCTCGCTCACTGACTCGCTACGCTCGGTCGTTCGACTGCGGCGAGCGGAAATGGCTTACGAACGGGGCGGAGATTTCCTGGAAGATGCCAGGAAGATACTTAACAGGGAAGTGAGAGGGCCGCGGCAAAGCCGTTTTTCCATAGGCTCCGCCCCCCTGACAAGCATCACGAAATCTGACGCTCAAATCAGTGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCTGGCGGCTCCCTCGTGCGCTCTCCTGTTCCTGCCTTTCGGTTTACCGGTGTCATTCCGCTGTTATGGCCGCGTTTGTCTCATTCCACGCCTGACACTCAGTTCCGGGTAGGCAGTTCGCTCCAAGCTGGACTGTATGCACGAACCCCCCGTTCAGTCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGAAAGACATGCAAAAGCACCACTGGCAGCAGCCACTGGTAATTGATTTAGAGGAGTTAGTCTTGAAGTCATGCGCCGGTTAAGGCTAAACTGAAAGGACAAGTTTTGGTGACTGCGCTCCTCCAAGCCAGTTACCTCGGTTCAAAGAGTTGGTAGCTCAGAGAACCTTCGAAAAACCGCCCTGCAAGGCGGTTTTTTCGTTTTCAGAGCAAGAGATTACGCGCAGACCAAAACGATCTCAAGAAGATCATCTTATTAAGGGGTCTGACGCTCAGTGGAACGAAAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCTTAGAAAAACTCATCGAGCATCAAATGAAACTGCAATTTATTCATATCAGGATTATCAATACCATATTTTTGAAAAAGCCGTTTCTGTAATGAAGGAGAAAACTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGGTATCGGTCTGCGATTCCGACTCGTCCAACATCAATACAACCTATTAATTTCCCCTCGTCAAAAATAAGGTTATCAAGTGAGAAATCACCATGAGTGACGACTGAATCCGGTGAGAATGGCAAAAGCTTATGCATTTCTTTCCAGACTTGTTCAACAGGCCAGCCATTACGCTCGTCATCAAAATCACTCGCATCAACCAAACCGTTATTCATTCGTGATTGCGCCTGAGCGAGACGAAATACGCGATCGCTGTTAAAAGGACAATTACAAACAGGAATCGAATGCAACCGGCGCAGGAACACTGCCAGCGCATCAACAATATTTTCACCTGAATCAGGATATTCTTCTAATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGGTGAGTAACCATGCATCATCAGGAGTACGGATAAAATGCTTGATGGTCGGAAGAGGCATAAATTCCGTCAGCCAGTTTAGTCTGACCATCTCATCTGTAACATCATTGGCAACGCTACCTTTGCCATGTTTCAGAAACAACTCTGGCGCATCGGGCTTCCCATACAATCGATAGATTGTCGCACCTGATTGCCCGACATTATCGCGAGCCCATTTATACCCATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTCGAGCAAGACGTTTCCCGTTGAATATGGCTCATAACACCCCTTGTATTACTGTTTATGTAAGCAGACAGTTTTATTGTTCATGATGATATATTTTTATCTTGTGCAATGTACATCAGAGATTTTGAGACACAACCAATTATTGAAGGCCTCCCTAACGGGGGGCCTTTTTTTGTTTCTGGTCTCCC'
seq_formatted = seq[-prefix_len:] + seq[0:-prefix_len]

# Make the part list with the positions and orientations
part_list = []

# AmtR unit
part_list.append({ 'part_name' : 'PTac',
	               'dir'       : 'F',
	               'seq_idx'   : 19 + prefix_len })
part_list.append({ 'part_name' : 'PTet',
	               'dir'       : 'F',
	               'seq_idx'   : 90 + prefix_len })
part_list.append({ 'part_name' : 'BydvJ',
	               'dir'       : 'F',
	               'seq_idx'   : 153 + prefix_len })
part_list.append({ 'part_name' : 'AmtR',
	               'dir'       : 'F',
	               'seq_idx'   : 263 + prefix_len,
	               'char_K'    : 0.169953,
	               'char_n'    : 1.319126,
	               'char_max'  : 13.18696,
	               'char_min'  : 0.316394,
	               'T1_1_FPKM' : 0.0120263,
	               'T1_2_FPKM' : 0.0114438,
	               'T2_1_FPKM' : 0.046604,
	               'T2_2_FPKM' : 0.0552453,
	               'T3_1_FPKM' : 0.0966424,
	               'T3_2_FPKM' : 0.0683571,
	               'T4_1_FPKM' : 0.118099,
	               'T4_2_FPKM' : 0.119218,
	               'T5_1_FPKM' : 0.0330196,
	               'T5_2_FPKM' : 0.0426507,
	               'T6_1_FPKM' : 0.0558763,
	               'T6_2_FPKM' : 0.0457136,
	               'T7_1_FPKM' : 0.235475,
	               'T7_2_FPKM' : 0.0810696,
	               'T8_1_FPKM' : 0.341767,
	               'T8_2_FPKM' : 0.114558 })
part_list.append({ 'part_name' : 'L3S2P55',
	               'dir'       : 'F',
	               'seq_idx'   : 932 + prefix_len })

# LitR unit
part_list.append({ 'part_name' : 'PBAD',
	               'dir'       : 'F',
	               'seq_idx'   : 1020 + prefix_len })
part_list.append({ 'part_name' : 'PTet',
	               'dir'       : 'F',
	               'seq_idx'   : 1339 + prefix_len })
part_list.append({ 'part_name' : 'PlmJ',
	               'dir'       : 'F',
	               'seq_idx'   : 1402 + prefix_len })
part_list.append({ 'part_name' : 'LitR',
	               'dir'       : 'F',
	               'seq_idx'   : 1517 + prefix_len,
	               'char_K'    : 0.138481,
	               'char_n'    : 1.542546,
	               'char_max'  : 10.20081,
	               'char_min'  : 0.263379,
	               'T1_1_FPKM' : 0.00382445,
	               'T1_2_FPKM' : 0.00235393,
	               'T2_1_FPKM' : 0.0092055,
	               'T2_2_FPKM' : 0.00292538,
	               'T3_1_FPKM' : 0.404388,
	               'T3_2_FPKM' : 0.302069,
	               'T4_1_FPKM' : 0.296094,
	               'T4_2_FPKM' : 0.438252,
	               'T5_1_FPKM' : 1.14506,
	               'T5_2_FPKM' : 0.672097,
	               'T6_1_FPKM' : 0.153965,
	               'T6_2_FPKM' : 0.89729,
	               'T7_1_FPKM' : 1.21826,
	               'T7_2_FPKM' : 1.15901,
	               'T8_1_FPKM' : 1.04427,
	               'T8_2_FPKM' : 1.44773 })
part_list.append({ 'part_name' : 'L3S2P24',
	               'dir'       : 'F',
	               'seq_idx'   : 2123 + prefix_len })

# BM3R1 unit
part_list.append({ 'part_name' : 'PBAD',
	               'dir'       : 'F',
	               'seq_idx'   : 2209 + prefix_len })
part_list.append({ 'part_name' : 'SarJ',
	               'dir'       : 'F',
	               'seq_idx'   : 2517 + prefix_len })
part_list.append({ 'part_name' : 'BM3R1',
	               'dir'       : 'F',
	               'seq_idx'   : 2633 + prefix_len,
	               'char_K'    : 0.627955,
	               'char_n'    : 2.945651,
	               'char_max'  : 2.168494,
	               'char_min'  : 0.019056,
	               'T1_1_FPKM' : 0.00176498,
	               'T1_2_FPKM' : 0.000980959,
	               'T2_1_FPKM' : 0.00145981,
	               'T2_2_FPKM' : 0.000935148,
	               'T3_1_FPKM' : 0.00554415,
	               'T3_2_FPKM' : 0.00315491,
	               'T4_1_FPKM' : 0.00387313,
	               'T4_2_FPKM' : 0.00456769,
	               'T5_1_FPKM' : 0.142254,
	               'T5_2_FPKM' : 0.136984,
	               'T6_1_FPKM' : 0.0174949,
	               'T6_2_FPKM' : 0.114058,
	               'T7_1_FPKM' : 0.0237474,
	               'T7_2_FPKM' : 0.145217,
	               'T8_1_FPKM' : 0.0216352,
	               'T8_2_FPKM' : 0.163363 })
part_list.append({ 'part_name' : 'L3S2P11',
	               'dir'       : 'F',
	               'seq_idx'   : 3212 + prefix_len })

# SrpR unit
part_list.append({ 'part_name' : 'PBM3R1',
	               'dir'       : 'F',
	               'seq_idx'   : 3288 + prefix_len })
part_list.append({ 'part_name' : 'PAmtR',
	               'dir'       : 'F',
	               'seq_idx'   : 3354 + prefix_len })
part_list.append({ 'part_name' : 'RiboJ10',
	               'dir'       : 'F',
	               'seq_idx'   : 3406 + prefix_len })
part_list.append({ 'part_name' : 'SrpR',
	               'dir'       : 'F',
	               'seq_idx'   : 3515 + prefix_len,
	               'char_K'    : 0.43338,
	               'char_n'    : 2.926147,
	               'char_max'  : 5.538373,
	               'char_min'  : 0.011192,
	               'T1_1_FPKM' : 0.281223,
	               'T1_2_FPKM' : 0.548636,
	               'T2_1_FPKM' : 0.0613146,
	               'T2_2_FPKM' : 0.0480011,
	               'T3_1_FPKM' : 0.0633603,
	               'T3_2_FPKM' : 0.0620076,
	               'T4_1_FPKM' : 0.0504884,
	               'T4_2_FPKM' : 0.0561277,
	               'T5_1_FPKM' : 0.242214,
	               'T5_2_FPKM' : 0.332731,
	               'T6_1_FPKM' : 0.0458752,
	               'T6_2_FPKM' : 0.0567179,
	               'T7_1_FPKM' : 0.032508,
	               'T7_2_FPKM' : 0.0739426,
	               'T8_1_FPKM' : 0.0432432,
	               'T8_2_FPKM' : 0.0557725 })
part_list.append({ 'part_name' : 'ECK120029600',
	               'dir'       : 'F',
	               'seq_idx'   : 4157 + prefix_len })

# PhlF unit
part_list.append({ 'part_name' : 'PSrpR',
	               'dir'       : 'F',
	               'seq_idx'   : 4266 + prefix_len })
part_list.append({ 'part_name' : 'PLitR',
	               'dir'       : 'F',
	               'seq_idx'   : 4348 + prefix_len })
part_list.append({ 'part_name' : 'RiboJ53',
	               'dir'       : 'F',
	               'seq_idx'   : 4400 + prefix_len })
part_list.append({ 'part_name' : 'PhlF',
	               'dir'       : 'F',
	               'seq_idx'   : 4512 + prefix_len,
	               'char_K'    : 0.560868,
	               'char_n'    : 3.876935,
	               'char_max'  : 17.22043,
	               'char_min'  : 0.073805,
	               'T1_1_FPKM' : 0.0826387,
	               'T1_2_FPKM' : 0.0860746,
	               'T2_1_FPKM' : 0.0629795,
	               'T2_2_FPKM' : 0.088388,
	               'T3_1_FPKM' : 0.0194236,
	               'T3_2_FPKM' : 0.0214471,
	               'T4_1_FPKM' : 0.0157426,
	               'T4_2_FPKM' : 0.0202562,
	               'T5_1_FPKM' : 0.0114903,
	               'T5_2_FPKM' : 0.0188559,
	               'T6_1_FPKM' : 0.0878957,
	               'T6_2_FPKM' : 0.0270693,
	               'T7_1_FPKM' : 0.166864,
	               'T7_2_FPKM' : 0.0408277,
	               'T8_1_FPKM' : 0.184441,
	               'T8_2_FPKM' : 0.0510033 })
part_list.append({ 'part_name' : 'ECK120033737',
	               'dir'       : 'F',
	               'seq_idx'   : 5115 + prefix_len })

# YFP unit
part_list.append({ 'part_name' : 'PPhlF',
	               'dir'       : 'F',
	               'seq_idx'   : 5191 + prefix_len })
part_list.append({ 'part_name' : 'RiboJ',
	               'dir'       : 'F',
	               'seq_idx'   : 5242 + prefix_len })
part_list.append({ 'part_name' : 'YFP',
	               'dir'       : 'F',
	               'seq_idx'   : 5343 + prefix_len,
	               'T1_1_FPKM' : 0.00508732,
	               'T1_2_FPKM' : 0.0053802,
	               'T2_1_FPKM' : 0.00249333,
	               'T2_2_FPKM' : 0.00384602,
	               'T3_1_FPKM' : 0.608534,
	               'T3_2_FPKM' : 0.650854,
	               'T4_1_FPKM' : 0.726605,
	               'T4_2_FPKM' : 0.747765,
	               'T5_1_FPKM' : 0.848437,
	               'T5_2_FPKM' : 1.18982,
	               'T6_1_FPKM' : 0.0355789,
	               'T6_2_FPKM' : 0.626831,
	               'T7_1_FPKM' : 0.243758,
	               'T7_2_FPKM' : 1.25245,
	               'T8_1_FPKM' : 0.205018,
	               'T8_2_FPKM' : 0.853487,
	               'T1_1_REU'  : 0.024889,
	               'T2_1_REU'  : 0.018281,
	               'T3_1_REU'  : 11.56606,
	               'T4_1_REU'  : 11.39615,
	               'T5_1_REU'  : 13.86618,
	               'T6_1_REU'  : 0.405305,
	               'T7_1_REU'  : 0.389572,
	               'T8_1_REU'  : 0.317202,
	               'T1_2_REU'  : 0.020169,
	               'T2_2_REU'  : 0.019225,
	               'T3_2_REU'  : 8.287373,
	               'T4_2_REU'  : 7.503886,
	               'T5_2_REU'  : 20.58717,
	               'T6_2_REU'  : 6.654322,
	               'T7_2_REU'  : 6.364841,
	               'T8_2_REU'  : 6.043894 })
part_list.append({ 'part_name' : 'L3S2P21',
	               'dir'       : 'F',
	               'seq_idx'   : 6063 + prefix_len })

# Make the circuit variant
circuit_lib.new_variant('1', part_list, seq=seq_formatted, attribs=None)

# Save to file
circuit_lib.save('alec_circuit_lib.txt')

######################################################################
# PROCESS AND CLEAN TRACE FILES
######################################################################

total_reads = {'1_1' : 2090122.0,
               '1_2' : 3764262.0,
               '2_1' : 2678407.0,
               '2_2' : 11774747.0,
               '3_1' : 2185618.0,
               '3_2' : 3178283.0,
               '4_1' : 9040281.0,
               '4_2' : 5147585.0,
               '5_1' : 1339545.0,
               '5_2' : 2530869.0,
               '6_1' : 2413638.0,
               '6_2' : 20048541.0,
               '7_1' : 100338.0,
               '7_2' : 3232192.0,
               '8_1' : 87257.0,
               '8_2' : 8187654.0}

def load_trace (trace_filename, prefix_len):
	# WARNING: no error checking - data must be in order
	trace = []
	found_data = False
	with open(trace_filename, 'r') as f:
		for line in f:
			if found_data == False:
				# Check to see if reached data
				if line[0:3] == 'Ref':
					found_data = True
			else:
				if len(line.split('\t')) == 2: 
					# Add the data to the trace if there is some
					trace.append(float(line.split('\t')[1]))
	# Return the trace with prefix length concatenated to start
	return trace[-prefix_len:] + trace[0:-prefix_len]

def clean_trace_file (fwd_filename, rev_filename, output_filename, prefix_len, tot_reads, norm_output_filename):
	fwd_trace = load_trace(fwd_filename, prefix_len)
	rev_trace = load_trace(rev_filename, prefix_len)
	# Save the traces with (bp, +, -) columns
	f_out = open(output_filename, 'w')
	f_out_norm = open(norm_output_filename, 'w')
	f_out.write('Pos,+,-\n')
	f_out_norm.write('Pos,+,-\n')
	for idx in range(len(fwd_trace)):
		f_out.write(str(idx) + ',' + str(fwd_trace[idx]) + ',' + str(rev_trace[idx]) + '\n')
		f_out_norm.write(str(idx) + ',' + str(float(fwd_trace[idx])/tot_reads) + ',' + str(float(rev_trace[idx])/tot_reads) + '\n')
	f_out.close()
	f_out_norm.close()

clean_trace_file ('./_raw/Circuit_0x58v50_T1_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T1_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T1_1.txt', prefix_len, total_reads['1_1'],
	              './_clean/alec_circuit_physical_trace_norm_T1_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T1_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T1_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T1_2.txt', prefix_len, total_reads['1_2'],
	              './_clean/alec_circuit_physical_trace_norm_T1_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T2_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T2_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T2_1.txt', prefix_len, total_reads['2_1'],
	              './_clean/alec_circuit_physical_trace_norm_T2_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T2_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T2_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T2_2.txt', prefix_len, total_reads['2_2'],
	              './_clean/alec_circuit_physical_trace_norm_T2_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T3_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T3_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T3_1.txt', prefix_len, total_reads['3_1'],
	              './_clean/alec_circuit_physical_trace_norm_T3_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T3_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T3_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T3_2.txt', prefix_len, total_reads['3_2'],
	              './_clean/alec_circuit_physical_trace_norm_T3_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T4_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T4_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T4_1.txt', prefix_len, total_reads['4_1'],
	              './_clean/alec_circuit_physical_trace_norm_T4_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T4_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T4_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T4_2.txt', prefix_len, total_reads['4_2'],
	              './_clean/alec_circuit_physical_trace_norm_T4_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T5_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T5_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T5_1.txt', prefix_len, total_reads['5_1'],
	              './_clean/alec_circuit_physical_trace_norm_T5_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T5_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T5_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T5_2.txt', prefix_len, total_reads['5_2'],
	              './_clean/alec_circuit_physical_trace_norm_T5_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T6_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T6_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T6_1.txt', prefix_len, total_reads['6_1'],
	              './_clean/alec_circuit_physical_trace_norm_T6_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T6_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T6_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T6_2.txt', prefix_len, total_reads['6_2'],
	              './_clean/alec_circuit_physical_trace_norm_T6_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T7_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T7_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T7_1.txt', prefix_len, total_reads['7_1'],
	              './_clean/alec_circuit_physical_trace_norm_T7_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T7_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T7_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T7_2.txt', prefix_len, total_reads['7_2'],
	              './_clean/alec_circuit_physical_trace_norm_T7_2.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T8_1.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T8_1.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T8_1.txt', prefix_len, total_reads['8_1'],
	              './_clean/alec_circuit_physical_trace_norm_T8_1.txt')

clean_trace_file ('./_raw/Circuit_0x58v50_T8_2.cleaned.merged.duplicates_marked.forward_strand.physical_coverage.metrics', 
	              './_raw/Circuit_0x58v50_T8_2.cleaned.merged.duplicates_marked.reverse_strand.physical_coverage.metrics',
	              './_clean/alec_circuit_physical_trace_T8_2.txt', prefix_len, total_reads['8_2'],
	              './_clean/alec_circuit_physical_trace_norm_T8_2.txt')

######################################################################
# PROCESS AND CLEAN TRANSCRIPTOMIC DATA
######################################################################

def load_tx_file (tx_filename):
	tx_data = {}
	# Open file and ignore headers
	tx_reader = csv.reader(open(tx_filename, 'rb'), delimiter='\t')
	header1 = next(tx_reader, None)
	header2 = next(tx_reader, None)
	# Load the data
	for row in tx_reader:
		tx_data[row[0]] = float(row[1])
	return tx_data

def compress_tx_dict (tx_dict_list):
	# Given a list of dictionaries of the tx data will compress into single matrix
	keys = tx_dict_list[0].keys()
	tx_matrix = np.zeros((len(keys),len(tx_dict_list)))
	i = 0
	for cur_key in keys:
		tx_matrix[i,:] = [x[cur_key] for x in tx_dict_list]
		i += 1
	return keys, tx_matrix

def save_tx_matrix (tx_keys, tx_matrix, output_filename):
	f_out = open(output_filename, 'w')
	f_out.write('gene')
	for i in range(np.size(tx_matrix, 1)):
		f_out.write(',state_'+ str(i+1))
	f_out.write('\n')
	for k in range(np.size(tx_matrix, 0)):
		f_out.write(tx_keys[k])
		for i in range(np.size(tx_matrix, 1)):
			f_out.write(','+ str(tx_matrix[k,i]))
		f_out.write('\n')
	f_out.close()

tx_data = []
tx_data.append(load_tx_file('./_raw/Circuit_sample_T1_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T2_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T3_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T4_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T5_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T6_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T7_1.FPKM'))
tx_data.append(load_tx_file('./_raw/Circuit_sample_T8_1.FPKM'))
tx_keys, tx_matrix = compress_tx_dict(tx_data)
save_tx_matrix(tx_keys, tx_matrix, './_clean/alec_circuit_tx_1.csv')

tx_bad_data = []
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T1_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T2_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T3_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T4_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T5_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T6_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T7_2.FPKM'))
tx_bad_data.append(load_tx_file('./_raw/Circuit_sample_T8_2.FPKM'))
tx_bad_keys, tx_bad_matrix = compress_tx_dict(tx_bad_data)
save_tx_matrix(tx_bad_keys, tx_bad_matrix, './_clean/alec_circuit_tx_2.csv')
