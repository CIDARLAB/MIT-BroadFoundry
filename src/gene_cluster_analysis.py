#!/usr/bin/env python
"""
Gene Cluster Library Analysis Tools
===================================

    Gene Cluster Analysis Tools are a set of methods to robustly perform 
    general queries of gene cluster libraries.
"""
#    Gene Cluster Tools for Python
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import sys
if sys.version_info[:2] < (2, 6):
    m = "Python version 2.6 or later is required for Gene Cluster Tools (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

import csv
import gene_cluster_library as gcl

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def transcriptional_units (gcl):
	"""Extract all valid transcriptional units from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    Returns
    -------
    units: dict(list([PromoterID, TerminatorID]))
        All valid transcriptional units i.e., Promoter -> Terminator. Dictionary 
        indexed by variant name with elements stored as a list of promoter, 
        terminator pairs (also stored as a list).
	"""
	# Find all promoters in the library
	p_insts = gcl.find_part_type_instances('Promoter')
	# Find all the downstream terminators from the promoters
	t_insts = gcl.find_next_part_idxs(p_insts, part_type='Terminator')
	# Combine the data to give the full transcriptional unit boundaries
	units = {}
	for v_key in p_insts.keys():
		if v_key in t_insts.keys():
			cur_p_data = p_insts[v_key]
			cur_t_data = t_insts[v_key]
			for idx in range(len(cur_p_data)):
				if cur_t_data[idx] != None:
					if v_key not in units.keys():
						units[v_key] = []
					units[v_key].append([cur_p_data[idx], cur_t_data[idx]])
	return units

def monocistronic_units (gcl):
	"""Extract all monocistronic transcriptional units from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    Returns
    -------
    units: dict(list([PromoterID, TerminatorID]))
        All monocistronic transcriptional units. Dictionary indexed by variant name 
        with elements stored as a list of promoter, terminator pairs (also stored as 
        a list).
	"""
	return polycistronic_units(gcl, number_of_CDSs=1)

def polycistronic_units (gcl, number_of_CDSs=None):
	"""Extract all polycistronic transcriptional units from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    number_of_CDSs : int (default=None)
    	Number of CDS parts required between promoter and terminator to be included.
    	If == None then all counts > 1 are included.

    Returns
    -------
    units: dict(list([PromoterID, TerminatorID]))
        All polycistronic transcriptional units. Dictionary indexed by variant name 
        with elements stored as a list of promoter, terminator pairs (also stored as 
        a list).
	"""
	# Find all promoters in the library
	p_insts = gcl.find_part_type_instances('Promoter')
	# Find all the downstream terminators from the promoters
	t_insts = gcl.find_next_part_idxs(p_insts, part_type='Terminator')
	# For each valid pair, count number of CDS part types between only include if == 1
	units = {}
	for v_key in p_insts.keys():
		if v_key in t_insts.keys():
			cur_p_data = p_insts[v_key]
			cur_t_data = t_insts[v_key]
			for idx in range(len(cur_p_data)):
				cur_p = cur_p_data[idx]
				cur_t = cur_t_data[idx]
				# Cycle through all parts between promoter and terminator and count CDSs
				cds_count = 0
				for cur_part_idx in range(cur_p+1,cur_t):
					cur_part_name = gcl.variants[v_key]['part_list'][cur_part_idx]
					if gcl.parts[cur_part_name]['type'] == 'CDS':
						cds_count += 1
				if number_of_CDSs == None:
					if cds_count > 1:
						if v_key not in units.keys():
							units[v_key] = []
						units[v_key].append([cur_p, cut_t])
				elif number_of_CDSs == cds_count:
					if v_key not in units.keys():
							units[v_key] = []
						units[v_key].append([cur_p, cut_t])
	return units

def extract_ranges_for_transcriptional_units (gcl, units):
	"""Extract sequences ranges for a set of transcriptional units.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    units : dict(list([start index, end index]))
    	Transcriptional unit start and end locations to extract.

    Returns
    -------
    ranges: dict(list([start seq idx, end seq index]))
        Sequence ranges for each transcriptional unit in the library. Stored as 
        a dictionary keyed on variant and then a list of tuples containing the 
        start and end sequence locations of the transcriptional units.
	"""
	# Cycle through all promoter/terminator pairs and find start/end indexes
	ranges = {}
	for v_key in units.keys():
		for el in units[v_key]:
			p_start_idx = gcl.variants[v_key]['part_list'][el[0]]['seq_idx']
			p_end_idx = p_start_idx + gcl.variants[v_key]['part_list'][el[0]]['seq_len']
			t_start_idx = gcl.variants[v_key]['part_list'][el[0]]['seq_idx']
			t_end_idx = t_start_idx + gcl.variants[v_key]['part_list'][el[0]]['seq_len']
			if v_key not in ranges.keys():
				ranges[v_key] = []
			# We return the full length of the transcriptional unit promoter start -> terminator end
			ranges[v_key].append([p_start_idx, t_end_idx])
	return ranges

def divergent_promoters (gcl):
	"""Extract all divergent promoters from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    Returns
    -------
    locations: dict(list(list))
        Locations of divergent promoters in the library. Stored as dictionary keyed
        on variant and then a list of tuples containing the locations of the 
        promoter pairs.
	"""
	return None

def convergent_promoters (gcl):
	"""Extract all convergent promoters from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    Returns
    -------
    locations: dict(list(list))
        Locations of convergent promoters in the library. Stored as dictionary keyed
        on variant and then a list of tuples containing the locations of the 
        promoter pairs.
	"""
	return None

def double_promoters (gcl):
	"""Extract all double promoters from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    Returns
    -------
    locations: dict(list(list))
        Locations of double promoters in the library. Stored as dictionary keyed
        on variant and then a list of tuples containing the locations of the 
        promoters.
	"""
	return None

def double_terminators (gcl):
	"""Extract all double terminators from a GeneClusterLibrary.

    Parameters
    ----------
    gcl : GeneClusterLibrary
        Gene cluster library to consider.

    Returns
    -------
    locations: dict(list(list))
        Locations of double terminators in the library. Stored as dictionary keyed
        on variant and then a list of tuples containing the locations of the 
        terminators.
	"""
	return None
