#!/usr/bin/env python
"""
Gene Cluster Library
====================

    Gene Cluster Library is a small library designed for creating, querying
    and saving information related to refactored gene clusters. It provides
    data structures and functions to manipulate light-weight text based
    files that contain all sequence, part and cluster variant information,
    in addition to search capabilities to facilitate easy querying of data
    during analysis. It can output clusters to PigeonCAD format to generate
    a visual representation or matplotlib from vector based figures.
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
import datetime

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

###############################################################################
# MAIN OBJECT FOR GENE CLUSTER LIBRARIES
###############################################################################

class GeneClusterLibrary:
	"""Class encapsulating the concept of a gene cluster and potential variants.
	"""

	# Static class variable holding standard mapping of part type to PigeonCAD type.
	STD_PIGEON_MAP = {'Promoter':'p',
	                  'CDS':'c',
	                  'RBS':'r',
	                  'Terminator':'t',
	                  'Scar':'x',
	                  'Spacer':'s'}

	STD_GENBANK_MAP = {'Promoter':'promoter',
	                  'CDS':'CDS',
	                  'RBS':'RBS',
	                  'Terminator':'terminator',
	                  'Scar':'misc_feature',
	                  'Spacer':'misc_feature',
	                  'Ribozyme':'misc_feature'}

	# Static class variable listing the standard parts generally used.
	STD_PART_TYPES = ['Promoter',
	                  'RBS',
	                  'CDS',
	                  'Terminator',
	                  'Scar',
	                  'Spacer',
	                  'Ribozyme']

	def __init__(self):
		"""Constructor to generate an empty GeneClusterLibrary.
		"""
		self.parts = {}
		self.variants = {}

	def clear (self):
		"""Clear all data associated with the GeneClusterLibrary.
		"""
		self.parts = {}
		self.variants = {}

	def new_part (self, name, part_type, seq, attribs=None):
		"""Create a new part to use in the variants.

		Parts are the pieces that make up the gene clusters. A part consists
		of a name, type and sequence, and can contain additional attributes associated
		with it e.g., performance information such as REU. All attributes are
		held as a dict to enable easy and fast keyed access. Parts that have
		been defined can be accessed using the 'parts' instance variable.

	    Parameters
	    ----------
	    name : string
	        Part name.

	    part_type : string
	    	Part type (see STD_PART_TYPES for examples)

	    seq : string
	    	Part sequence.

	    attribs : dict (default=None)
	    	Part attributes

	    Returns
	    -------
	    the_part: dict
	        The new part
		"""
		self.parts[name] = {}
		the_part = self.parts[name]
		the_part['seq'] = seq
		the_part['type'] = part_type
		if attribs != None:
			# Handle potential changes in part attributes
			for k in attribs.keys():
				the_part[k] = attribs[k]
		return the_part

	def new_variant (self, name, part_list, seq=None, attribs=None):
		"""Create a new gene cluster variant.

		Variants are different potential instances of the same gene 
		cluster by using different sets of parts or a different genetic
		architeture. To create a variant it must be given a name and
		a list of parts in the order required. Each element in the part
		list must itself be a dict containing the part name ('part_name')
		and direction ('dir') keys. If only this information is given then
		a sequence is automatically generated from the part seqeunces held
		in the gene cluster. Alternatively, if there are other sequences
		in the construct that need to be included, the actual sequence
		can be given and then the part list elements must also contain the
		'seq_idx' that contains the index at which the start of the element
		is in the provided sequence.

	    Parameters
	    ----------
	    name : string
	        Variant name.

	    part_list : list(dict)
	    	Part list that defines the variant. The list contains dicts
	    	which define the specific ordering of the elements. Each dict
	    	should contain the following keys:
	    		'part_name' : the name of the part.
	    		'dir' : direction of part either forward ('F') or reverse ('R').
	    		'seq_idx' : start index in the sequence (optional and not
	    			        needed if sequence automatically generated). 

	    seq : string (default=None)
	    	Sequence of the variant if not fully annotated by the parts.
	    	If None is provided the sequence will be built automatically
	    	from the available part data.

	    attribs : dict (default=None)
	    	Variant attributes.

	    Returns
	    -------
	    the_part: dict
	        The new variant
		"""
		the_var = {}
		if seq == None:
			# Part list will not include seq_idxs (build them here)
			the_seq = ''
			the_var['part_list'] = []
			for el in part_list:
				part_list_el = {}
				part_list_el['part_name'] = el['part_name']
				part_list_el['dir'] = el['dir']
				part_list_el['seq_idx'] = len(the_seq)
				part_list_el['seq_len'] = len(self.parts[el['part_name']]['seq'])
				# Update the sequence
				the_seq = the_seq + self.parts[el['part_name']]['seq']
				# Add additional variant-part-level attributes
				for k in list(set(el.keys()) - set(['part_name', 'dir', 'seq_idx'])):
					part_list_el[k] = el[k]
				# Add part to list
				the_var['part_list'].append(part_list_el)
			the_var['seq'] = the_seq
		else:
			the_var['part_list'] = []
			for el in part_list:
				part_list_el = {}
				part_list_el['part_name'] = el['part_name']
				part_list_el['dir'] = el['dir']
				part_list_el['seq_idx'] = el['seq_idx']
				part_list_el['seq_len'] = len(self.parts[el['part_name']]['seq'])
				# Add additional variant-part-level attributes
				for k in list(set(el.keys()) - set(['part_name', 'dir', 'seq_idx'])):
					part_list_el[k] = el[k]
				# Add part to list
				the_var['part_list'].append(part_list_el)
			the_var['seq'] = seq
		if attribs != None:
			# Add additional variant-level attributes
			for k in attribs.keys():
				the_var[k] = attribs[k]
		self.variants[name] = the_var
		return the_var

	def variant_data (self, variant):
		"""Get method for variant data from a name

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    Returns
	    -------
	    vairant_data: dict
	        Dictionary of the variant data.
		"""
		return self.variants[variant]

	def variant_part_list (self, variant):
		"""Get method for variant part list from a name

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    Returns
	    -------
	    part_list: list
	        Part list for the variant.
		"""
		return self.variants[variant]['part_list']

	def variant_part_idx_name (self, variant, part_idx):
		"""Get method for part name from variant and part index

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    Returns
	    -------
	    part_name: string
	        Name of the part.
		"""
		return self.variants[variant]['part_list'][part_idx]['part_name']

	def variant_part_idx_type (self, variant, part_idx):
		"""Get method for part type from variant and part index

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    Returns
	    -------
	    part_type: string
	        Type of the part.
		"""
		return self.parts[self.get_variant_part_idx_name(variant, part_idx)]['type']

	def variant_part_idx_dir (self, variant, part_idx):
		"""Get method for part direction from variant and part index

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    Returns
	    -------
	    dir: string
	        Direction of part ('F' = Forward, 'R' = Reverse).
		"""
		return self.variants[variant]['part_list'][part_idx]['dir']

	def variant_part_idx_seq_idx (self, variant, part_idx):
		"""Get method for start bp from variant and part index

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    Returns
	    -------
	    seq_idx: int
	        Index of the start bp on +ve strand of the part.
		"""
		return self.variants[variant]['part_list'][part_idx]['seq_idx']

	def variant_part_idx_seq_len (self, variant, part_idx):
		"""Get method for part sequence length from variant and part index

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    Returns
	    -------
	    seq_len: int
	        Length of the part.
		"""
		return self.variants[variant]['part_list'][part_idx]['seq_len']

	def variant_part_idx_attrib (self, variant, part_idx, attrib_key):
		"""Get method for a part's attribute from variant, part index and key

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    attrib_key : string
	    	Attribute key.

	    Returns
	    -------
	    attrib_value: undefined object
	        Attribute data.
		"""
		if attrib_key in self.variants[variant]['part_list'][part_idx].keys():
			return self.variants[variant]['part_list'][part_idx][attrib_key]
		else:
			return None

	def set_variant_part_idx_attrib (self, variant, part_idx, attrib_key, attrib_value):
		"""Set method for a part's attribute from variant, part index and key

	    Parameters
	    ----------
	    variant : string
	        Variant name.

	    part_idx : int
	    	Part index.

	    attrib_key : string
	    	Attribute key.

	    attrib_value : string
	    	Attribute value.
		"""
		self.variants[variant]['part_list'][part_idx][attrib_key] = attrib_value

	def part_data (self, part_name):
		"""Get method for variant data from a name

	    Parameters
	    ----------
	    part_name : string
	        Part name.

	    Returns
	    -------
	    part_data: dict
	        Dictionary of the part data.
		"""
		return self.parts[part_name]

	def part_attrib (self, part_name, attrib_key):
		"""Get method for a part's attribute from a part name and key

	    Parameters
	    ----------
	    part_name : string
	    	Part name.

	    attrib_key : string
	    	Attribute key.

	    Returns
	    -------
	    attrib_value: undefined object
	        Attribute data.
		"""
		if part_name in self.parts.keys():
			if attrib_key in self.parts[part_name].keys():
				return self.parts[part_name][attrib_key]
		return None

	def part_type (self, part_name):
		"""Get method for a part's type from a part name

	    Parameters
	    ----------
	    part_name : string
	    	Part name.

	    Returns
	    -------
	    part_type: string
	        Part type.
		"""
		if part_name in self.parts.keys():
			return self.parts[part_name]['type']
		return None

	def part_seq (self, part_name):
		"""Get method for a part's sequence from a part name

	    Parameters
	    ----------
	    part_name : string
	    	Part name.

	    Returns
	    -------
	    part_seq: string
	        Part sequence.
		"""
		if part_name in self.parts.keys():
			return self.parts[part_name]['seq']
		return None

	def variant_part_idx_start_bp (self, variant, part_idx):
		"""Find the start bp of a variant part (consdiers direction)

	    Parameters
	    ----------
	    variant : string
	    	Variant name.

	    part_name : string
	    	Part name.

	    Returns
	    -------
	    start_bp: int
	        Index (bp) of the start of the part.
		"""
		if part_idx == None:
			return None
		else:
			part_data = self.variants[variant]['part_list'][part_idx]
			if part_data['dir'] == 'F':
				return part_data['seq_idx']
			else:
				return part_data['seq_idx'] + part_data['seq_len'] 

	def variant_part_idx_end_bp (self, variant, part_idx):
		"""Find the end bp of a variant part (consdiers direction)

	    Parameters
	    ----------
	    variant : string
	    	Variant name.

	    part_name : string
	    	Part name.

	    Returns
	    -------
	    end_bp: int
	        Index (bp) of the end of the part.
		"""
		if part_idx == None:
			return None
		else:
			part_data = self.variants[variant]['part_list'][part_idx]
			if part_data['dir'] == 'F':
				return part_data['seq_idx'] + part_data['seq_len']
			else:
				return part_data['seq_idx'] 

	def find_part_instances (self, part_name):
		"""Find instances of a part in all the variants.

	    Parameters
	    ----------
	    part_name : string
	        Part to search for.

	    Returns
	    -------
	    found : dict(list)
	        A dictionary keyed by the variant name that gives a list of all indexes
	        at which the specific part was found.
		"""
		found = {}
		for v in self.variants.keys():
			idxs = []
			for p_idx in range(len(self.variants[v]['part_list'])):
				if self.variants[v]['part_list'][p_idx]['part_name'] == part_name:
					idxs.append(p_idx)
			if idxs != []:
				found[v] = idxs
		return found

	def find_part_type_instances (self, part_type):
		"""Find instances of a part in all the variants.

	    Parameters
	    ----------
	    part_type : string
	        Part type to search for.

	    Returns
	    -------
	    found : dict(list)
	        A dictionary keyed by the variant name that gives a list of all indexes
	        at which the part type was found.
		"""
		found = {}
		for v in self.variants.keys():
			idxs = []
			for p_idx in range(len(self.variants[v]['part_list'])):
				if self.parts[self.variants[v]['part_list'][p_idx]['part_name']]['type'] == part_type:
					idxs.append(p_idx)
			if idxs != []:
				found[v] = idxs
		return found

	def extract_seq_range (self, variant_name, part_idx, start_offset, end_offset):
		"""Extract sequence range from the start of a part in a specific variant.

	    Parameters
	    ----------
	    variant_name : string
	        Variant to use.

	    part_idx : string
	        Index of the part in the variant part list to use.

	    start_offset : int
	        Number of base pairs in the forward direction to include.

	    end_offset : int
	        Number of base pairs in the reverse direction to include.

	    Returns
	    -------
	    seq : string
	        The extracted sequence.
		"""
		the_part = self.variants[variant_name]['part_list'][part_idx]
		the_seq = self.variants[variant_name]['seq']
		the_seq_idx = the_part['seq_idx']
		the_seq_len = the_part['seq_len']
		if the_part['dir'] == 'F':
			start_idx = the_seq_idx-start_offset
			end_idx = the_seq_idx+end_offset
			return the_seq[start_idx:end_idx]
		else:
			start_idx = the_seq_idx+the_seq_len+start_offset
			end_idx = the_seq_idx+the_seq_len-end_offset
			return the_seq[end_idx:start_idx]

	def extract_seq_ranges (self, variant_insts, start_offset, end_offset):
		"""Extract sequence ranges from the start of a set of part instances.
		
	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant to use.

	    start_offset : int
	        Number of base pairs in the forward direction to include.

	    end_offset : int
	        Number of base pairs in the reverse direction to include.

	    Returns
	    -------
	    seqs : dict(list(int))
	        The extracted sequence.
		"""
		part_seqs = {}
		for v in variant_insts.keys():
			v_seqs = []
			for i in variant_insts[v]:
				v_seqs.append(self.extract_seq_range(v, i, start_offset, end_offset))
			part_seqs[v] = v_seqs
		return part_seqs

	def extract_seq_range_around_part (self, variant_name, part_idx, fwd_len, rev_len):
		"""Extract sequences that extend forward and backwards from the 
		edges of a part.

	    Parameters
	    ----------
	    variant_name : string
	        Variant to use.

	    part_idx : string
	        Index of the part in the variant part list to use.

	    fwd_len : int
	        Number of base pairs in the forward direction from the part 
	        end to include.

	    rev_len : int
	        Number of base pairs in the reverse direction from the part
	        start to include.

	    Returns
	    -------
	    seqs : dict(list(string))
	        The extracted sequence.
		"""
		part_len = len(self.parts[self.variants[variant_name]['part_list'][part_idx]['part_name']]['seq'])
		return self.extract_seq_range(variant_name, part_idx, part_len+fwd_len, rev_len)

	def extract_seq_ranges_around_part (self, variant_insts, fwd_len, rev_len):
		"""Extract sequences that extend forward and backwards from the 
		edges of a set of part instances.

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant to use.

	    fwd_len : int
	        Number of base pairs in the forward direction from the part 
	        end to include.

	    rev_len : int
	        Number of base pairs in the reverse direction from the part
	        start to include.

	    Returns
	    -------
	    seqs : dict(list(string))
	        The extracted sequence.
		"""
		part_seqs = {}
		for v in variant_insts.keys():
			v_seqs = []
			for i in variant_insts[v]:
				v_seqs.append(self.extract_seq_range_around_part(v, i, fwd_len, rev_len))
			part_seqs[v] = v_seqs
		return part_seqs

	def find_next_part_idx (self, variant_name, part_idx, part_type=None, next_count=1, dir_check=False):
		"""Find the next part in the design from specific part in a variant.

	    Parameters
	    ----------
	    variant_name : string
	        Variant name to use.

	    part_idx : int
	        Part index in variant part list to use.

	    part_type : string (default=None)
	        The type of part to search for in the forward direction. If None is
	        given then the next part irrespective of type is returned.

	    next_count : int
	        The number of parts to skip.

	    remove_nones : bool
	    	Flag as to whether potential None values should be included where no 
	    	proceeding part is found.

	    Returns
	    -------
	    idx : int
	        The next part index (or None if none exist)
		"""
		the_part = self.variants[variant_name]['part_list'][part_idx]
		cur_next_count = 1
		if the_part['dir'] == 'F':
			if part_type == None:
				next_part_idx = part_idx+next_count
				if next_part_idx < len(self.variants[variant_name]['part_list']):
					return next_part_idx
				else:
					return None
			else:
				for i in range(part_idx+1, len(self.variants[variant_name]['part_list'])):
					if self.parts[self.variants[variant_name]['part_list'][i]['part_name']]['type'] == part_type:
						if dir_check == True:
							if self.variants[variant_name]['part_list'][i]['dir'] == 'F':
								if cur_next_count == next_count:
									return i
								else:
									cur_next_count += 1
						else:
							if cur_next_count == next_count:
								return i
							else:
								cur_next_count += 1
		else:
			if part_type == None:
				next_part_idx = part_idx-next_count
				if next_part_idx >= 0:
					return next_part_idx
				else:
					return None
			else:
				for i in range(part_idx-1, -1, -1):
					if self.parts[self.variants[variant_name]['part_list'][i]['part_name']]['type'] == part_type:
						if dir_check == True:
							if self.variants[variant_name]['part_list'][i]['dir'] == 'R':
								if cur_next_count == next_count:
									return i
								else:
									cur_next_count += 1
						else:
							if cur_next_count == next_count:
								return i
							else:
								cur_next_count += 1
		return None

	def find_next_part_idxs (self, variant_insts, part_type=None, next_count=1, remove_nones=False, dir_check=False):
		"""Find the next parts in the design from variant part instances.

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant part instances to use.

	    part_type : string (default=None)
	        The type of part to search for in the forward direction. If None is
	        given then the next part irrespective of type is returned.

	    next_count : int
	        The number of parts to skip.

	    remove_nones : bool
	    	Flag as to whether potential None values should be included where no 
	    	proceeding part is found.

	    Returns
	    -------
	    part_idxs : dict(list(int))
	        Variant part instances of the found parts
		"""
		part_idxs = {}
		for v in variant_insts.keys():
			v_idxs = []
			for i in variant_insts[v]:
				v_idxs.append(self.find_next_part_idx(v, i, part_type, next_count, dir_check))
			part_idxs[v] = v_idxs
		if remove_nones:
			return self.__remove_nones_variant_data(part_idxs)
		else:
			return part_idxs

	def find_prev_part_idx (self, variant_name, part_idx, part_type=None, next_count=1, dir_check=False):
		"""Find the previous part in the design from specific part in a variant.

	    Parameters
	    ----------
	    variant_name : string
	        Variant name to use.

	    part_idx : int
	        Part index in variant part list to use.

	    part_type : string (default=None)
	        The type of part to search for in the reverse direction. If None is
	        given then the next part irrespective of type is returned.

	    next_count : int
	        The number of parts to skip.

	    remove_nones : bool
	    	Flag as to whether potential None values should be included where no 
	    	preceeding part is found.

	    Returns
	    -------
	    idx : int
	        The previous part index (or None if none exist)
		"""
		the_part = self.variants[variant_name]['part_list'][part_idx]
		cur_next_count = 1
		if the_part['dir'] == 'F':
			if part_type == None:
				next_part_idx = part_idx-next_count
				if next_part_idx >= 0:
					return next_part_idx
				else:
					return None
			else:
				for i in range(part_idx-1, -1, -1):
					if self.parts[self.variants[variant_name]['part_list'][i]['part_name']]['type'] == part_type:
						if dir_check == True:
							if self.variants[variant_name]['part_list'][i]['dir'] == 'F':
								if cur_next_count == next_count:
									return i
								else:
									cur_next_count += 1
						else:
							if cur_next_count == next_count:
								return i
							else:
								cur_next_count += 1
		else:
			if part_type == None:
				next_part_idx = part_idx+next_count
				if next_part_idx < len(self.variants[variant_name]['part_list']):
					return next_part_idx
				else:
					return None
			else:
				for i in range(part_idx+1, len(self.variants[variant_name]['part_list'])):
					if self.parts[self.variants[variant_name]['part_list'][i]['part_name']]['type'] == part_type:
						if dir_check == True:
							if self.variants[variant_name]['part_list'][i]['dir'] == 'F':
								if cur_next_count == next_count:
									return i
								else:
									cur_next_count += 1
						else:
							if cur_next_count == next_count:
								return i
							else:
								cur_next_count += 1
		return None
	
	def find_prev_part_idxs (self, variant_insts, part_type=None, next_count=1, remove_nones=False, dir_check=False):
		"""Find the previous parts in the design from variant part instances.

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant part instances to use.

	    part_type : string (default=None)
	        The type of part to search for in the reverse direction. If None is
	        given then the next part irrespective of type is returned.

	    next_count : int
	        The number of parts to skip.

	    remove_nones : bool
	    	Flag as to whether potential None values should be included where no 
	    	proceeding part is found.

	    Returns
	    -------
	    part_idxs : dict(list(int))
	        Variant part instances of the found parts
		"""
		part_idxs = {}
		for v in variant_insts.keys():
			v_idxs = []
			for i in variant_insts[v]:
				v_idxs.append(self.find_prev_part_idx(v, i, part_type, next_count, dir_check))
			part_idxs[v] = v_idxs
		if remove_nones:
			return self.__remove_nones_variant_data(part_idxs)
		else:
			return part_idxs

	def find_seq_idx_range (self, variant_name, part_idx):
		"""Find the start and end index in the sequence for a given part in
		in a specific variant.

	    Parameters
	    ----------
	    variant_name : string
	        Variant name to use.

	    part_idx : int
	        Part index in variant part list to use.

	    Returns
	    -------
	    range : (int, int)
	        Start and end index. If second index is lower than first then element
	        is in reverse orientation.
		"""
		if part_idx == None:
			None
		else:
			the_part = self.variants[variant_name]['part_list'][part_idx]
			if the_part['dir'] == 'F':
				return (the_part['seq_idx'], the_part['seq_idx']+the_part['seq_len'])
			else:
				return (the_part['seq_idx']+the_part['seq_len'], the_part['seq_idx'])

	def find_seq_idx_ranges (self, variant_insts):
		"""Find the start and end index in the sequences for a given set of 
		variant part instances.

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant part instances to use.

	    Returns
	    -------
	    ranges : dict(list((int, int)))
	        Dictionary for each variant that returns a list of start and end indexes. 
	        If second index is lower than first then element is in reverse orientation.
		"""
		part_ranges = {}
		for v in variant_insts.keys():
			v_ranges = []
			for i in variant_insts[v]:
				v_ranges.append(self.find_seq_idx_range(v, i))
			part_ranges[v] = v_ranges
		return part_ranges

	def extract_part_attribs (self, variant_insts, attrib_name, collate=False):
		"""Extract an attribute from variant part instances

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant part instances to use.

	    attrib_name : string
	    	Attribute name to extract.

	   	collate : bool
	   		Flag whether to collate all part variant data into a single list.

	    Returns
	    -------
	    part_attribs : dict(list) or list
	        Depending on if the attributes found are collated either a dict of lists
	        (non-collated) or single list (collated) is returned.
		"""
		part_attribs = {}
		for v in variant_insts.keys():
			v_attribs = []
			for i in variant_insts[v]:
				if i == None:
					v_attribs.append(None)
				else:
					v_attribs.append(self.variants[v]['part_list'][i][attrib_name])
			part_attribs[v] = v_attribs
		if collate:
			return self.__collate_variant_data_list(part_attribs)
		else:
			return part_attribs

	def extract_variant_attribs (self, variant_list, attrib_name, collate=False):
		"""Extract a specific attributes from a list of variant names

	    Parameters
	    ----------
	    variant_list : list(string)
	        List of variants to extract the attribute for.

	    attrib_name : string
	    	Attribute name to extract.

	   	collate : bool
	   		Flag whether to collate all part variant data into a single list.

	    Returns
	    -------
	    variant_attribs : dict(list) or list
	        Depending on if the attributes found are collated either a dict of lists
	        (non-collated) or single list (collated) is returned.
		"""
		variant_attribs = {}
		for v in variant_list:
			variant_attribs[v] = self.variants[v][attrib_name]
		if collate:
			return self.__collate_variant_data_value(variant_attribs)
		else:
			return variant_attribs

	def find_next_parts (self, variant_insts, n_next=1, dir_check=False):
		"""Extract part names up to a given distance forward from a current part.

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant part instances to use.

	   	n_next : int
	   		Number of next parts to extract e.g., 2 would extract all parts
	   		a distance of 1 and 2 away from the original part index.

	    Returns
	    -------
	    next_parts : list(string)
	        List of part_names.
		"""
		next_parts = []
		for i in range(1, n_next+1):
			part_insts = self.find_next_part_idxs(variant_insts, next_count=i, remove_nones=True, dir_check=dir_check)
			next_parts = next_parts + self.extract_part_names_from_idxs(part_insts)
		return next_parts

	def find_prev_parts (self, variant_insts, n_prev=1, dir_check=False):
		"""Extract part names up to a given distance in reverse from a current part.

	    Parameters
	    ----------
	    variant_insts : dict(list(int))
	        Variant part instances to use.

	   	n_next : int
	   		Number of next parts to extract e.g., 2 would extract all parts
	   		a distance of 1 and 2 away from the original part index.

	    Returns
	    -------
	    prev_parts : list(string)
	        List of part_names.
		"""
		prev_parts = []
		for i in range(1, n_prev+1):
			part_insts = self.find_prev_part_idxs(variant_insts, next_count=i, remove_nones=True, dir_check=dir_check)
			prev_parts = prev_parts + self.extract_part_names_from_idxs(part_insts)
		return prev_parts

	def extract_part_names_from_idxs (self, variant_insts):
		"""Extract part names from a list of variant instances.

	    Parameters
	    ----------
	    variant_insts : dict(list)
	        List of part variant instances to collate part names for.

	    Returns
	    -------
	    part_list : list(string)
	    	List of all part names in the variants.	    
		"""
		part_list = []
		for v in variant_insts.keys():
			for i in variant_insts[v]:
				part_list.append(self.variants[v]['part_list'][i]['part_name'])
		return part_list

	def full_variants_part_list (self, variants_list=None):
		"""Return list of all part names in all variants or subset of variants.

	    Parameters
	    ----------
	    variants_list : list(string) (default=None)
	        List of variants to collate all part names for. If None then all 
	        variants analysed.

	    Returns
	    -------
	    part_list : list(string)
	    	List of all part names in the variants.	    
		"""
		part_list = []
		v_keys = []
		if variants_list == None:
			v_keys = variants.keys()
		else:
			v_keys = variants_list
		for v in v_keys:
			for i in variants[v]['part_list']:
				part_list.append(i['part_name'])
		return part_list

	def load (self, filename):
		"""Load gene cluster information from a file.

	    Parameters
	    ----------
	    filename : string
	        Name of the file to load.
		"""
		# Clear the current object data
		self.clear()
		# Load entire file into list
		f = open(filename, 'r')
		raw_lines = f.readlines()
		# Strip empty lines and comments (#)
		lines = []
		for l in raw_lines:
			if len(l) != 0 and l[0] != '#' and l[0].lstrip() != '':
				lines.append(l.rstrip('\r\n'))
		# Work through file and load parts and variants
		reading_parts = False
		reading_variants = False
		l_idx = 0
		while l_idx < len(lines):
			cur_line = lines[l_idx]
			# Check which region currently in
			if cur_line == '>> PARTS':
				reading_parts = True
				reading_variants = False
			if cur_line == '>> VARIANTS':
				reading_parts = False
				reading_variants = True
			# Based on the current region type read the parts/variants
			if cur_line[0:2] == '> ':
				if reading_parts:
					part_name = cur_line[2::]
					part_type = lines[l_idx+1][6::]
					part_seq = lines[l_idx+2][5::]
					part_attribs = {}
					part_attribs_txt = lines[l_idx+3][9::].split(',')
					if part_attribs_txt != ['']:
						for a in part_attribs_txt:
							key_val = a.split('=')
							part_attribs[key_val[0]] = self.__make_float_if_needed(key_val[1])
					self.new_part(part_name, part_type, part_seq, attribs=part_attribs)
					l_idx += 3
				if reading_variants:
					var_name = cur_line[2::]
					var_seq = lines[l_idx+1][5::]
					var_part_list_txt = lines[l_idx+2][11::].split('; ')
					var_part_list = []
					for p in var_part_list_txt:
						p_split = p.split(',')
						cur_part = {}
						cur_part['part_name'] = p_split[0]
						cur_part['dir'] = p_split[1]
						cur_part['seq_idx'] = int(p_split[2])
						cur_part['seq_len'] = int(p_split[3])
						cur_part_attribs_txt = p_split[4::]
						for a in cur_part_attribs_txt:
							key_val = a.split('=')
							cur_part[key_val[0]] = self.__make_float_if_needed(key_val[1])
						var_part_list.append(cur_part)
					var_attribs = {}
					var_attribs_txt = lines[l_idx+3][9::].split(',')
					if var_attribs_txt != ['']:
						for a in var_attribs_txt:
							key_val = a.split('=')
							var_attribs[key_val[0]] = self.__make_float_if_needed(key_val[1])
					self.new_variant(var_name, var_part_list, seq=var_seq, attribs=var_attribs)
					l_idx += 3
			# Move to next line
			l_idx += 1
		f.close()

	def save (self, filename):
		"""Save gene cluster information to a file.

	    Parameters
	    ----------
	    filename : string
	        Name of the file to save gene cluster to.
		"""
		f = open(filename, 'w')
		# Header with library versioning
		f.write('# Creator: Gene Cluster Tools (version ' + __version__ + ')\n\n')
		# Write part information
		f.write('>> PARTS\n\n')
		for p in sorted(self.parts.keys()):
			f.write('> ' + str(p) + '\n')
			f.write('type: ' + str(self.parts[p]['type']) + '\n')
			f.write('seq: ' + str(self.parts[p]['seq']) + '\n')
			f.write('attribs: ')
			attribs_keys = sorted(list(set(self.parts[p].keys()) - set(['type', 'seq'])))
			if attribs_keys != []:
				attribs_str = ''
				for k in attribs_keys:
					attribs_str = attribs_str + str(k) + '=' + str(self.parts[p][k]) + ','
				f.write(attribs_str[:-1])
			f.write('\n\n')
		# Write variant information
		f.write('>> VARIANTS\n\n')
		for v in sorted(self.variants.keys()):
			f.write('> ' + str(v) + '\n')
			f.write('seq: ' + str(self.variants[v]['seq']) + '\n')
			# Generate the part lists
			f.write('part_list: ')
			parts_list_str = ''
			for p in self.variants[v]['part_list']:
				part_str = ''
				part_str = part_str + str(p['part_name']) + ',' + str(p['dir']) + ',' + str(p['seq_idx']) + ',' + str(p['seq_len']) + ','
				# Append attributes
				attribs_keys = sorted(list(set(p.keys()) - set(['part_name', 'dir', 'seq_idx', 'seq_len'])))
				for a in attribs_keys:
					part_str = part_str + str(a) + '=' + str(p[a]) + ','
				parts_list_str = parts_list_str + part_str[:-1] + '; '
			f.write(parts_list_str[:-2] + '\n')
			# Write variant attributes
			f.write('attribs: ')
			attribs_str = ''
			attribs_keys = sorted(list(set(self.variants[v].keys()) - set(['seq', 'part_list'])))
			for a in attribs_keys:
				attribs_str = attribs_str + str(a) + '=' + str(self.variants[v][a]) + ','
			f.write(attribs_str[:-1] + '\n\n')
		f.close()

	def save_variant_genbank (self, variant, filename, methylated=1, dna_type='linear', mapping=STD_GENBANK_MAP):
		"""Save gene cluster variant to a genbank format file.

	    Parameters
	    ----------
	    variant : string
	    	Variant to save to file

	    filename : string
	        Name of the file to save.
		"""
		f = open(filename, 'w')
		variant_name = variant + (' '*(24-len(variant)))
		var_seq = self.variants[variant]['seq'].upper()
		seq_len = len(var_seq)
		today = datetime.date.today()
		date_str = today.strftime('%d-%b-%Y')
		f.write('LOCUS       ' + variant_name +  str(seq_len) + ' bp ds-DNA     ' + dna_type + '       ' + date_str + '\n')
		f.write('DEFINITION  .\nACCESSION   \nVERSION     \nSOURCE      .\n  ORGANISM  .\nCOMMENT     \nCOMMENT     \nCOMMENT     ApEinfo:methylated:' + str(methylated) + '\n')
		f.write('FEATURES             Location/Qualifiers\n')
		# Save all the features as correct annotation
		for el in self.variants[variant]['part_list']:
			f.write(' '*5)
			part_name = el['part_name']
			part_type = self.parts[part_name]['type']
			part_gb_type = 'misc_feature'
			if part_type in mapping:
				part_gb_type = mapping[part_type]
			part_start_bp = str(el['seq_idx']+1)
			part_end_bp = str(el['seq_idx']+el['seq_len'])
			f.write(part_gb_type + ' '*(16-len(part_gb_type)))
			f.write(part_start_bp + '..' + part_end_bp + '\n')
			f.write('                     /label=' + part_name + '\n')
		# Save the sequence
		f.write('ORIGIN')
		for bp in range(0,len(var_seq),10):
			if bp % 60 == 0:
				bp_num = bp + 1
				# New line and add number at the start of the line
				f.write('\n' + ' '*(9-len(str(bp_num))) + str(bp_num))
			start_bp = bp
			end_bp = start_bp + 10
			if end_bp > len(var_seq):
				end_bp = len(var_seq)
			f.write(' ' + var_seq[start_bp:end_bp])
		f.write('\n//\n')
		f.close()

	def save_PigeonCAD (self, filename_prefix, mapping=STD_PIGEON_MAP, col_mapping=None):
		"""Save gene cluster variants to a set of PigeonCAD files for visualisation.

	    Parameters
	    ----------
	    filename_prefix : string
	        Prefix to use for each of the output files. Files will have format:
	        'filename_prefix + variant_name + .txt'

	    mapping : dict(string)
	    	Mapping to give a valid PigeonCAD string for a particular part type.

	    col_mapping : dict(int)
	    	Mapping to give a colour (integer 0-14) for each part in the output file.
		"""
		# Save each variant to separate PigeonCAD file
		for v in self.variants.keys():
			self.save_variant_PigeonCAD(v, filename_prefix + str(v) + '.txt', mapping=mapping, col_mapping=col_mapping)

	def save_variant_PigeonCAD (self, variant_name, filename, mapping=STD_PIGEON_MAP, col_mapping=None):
		"""Save single variant to a PigeonCAD file for visualisation.

	    Parameters
	    ----------
		variant_name : string
			Name of the variant to output.

	    filename : string
	        Name of the file to output to.

	    mapping : dict(string)
	    	Mapping to give a valid PigeonCAD string for a particular part type.

	    col_mapping : dict(int)
	    	Mapping to give a colour (integer 0-14) for each part in the output file.
		"""
		f_out = open(filename, 'w')
		part_list = self.variants[variant_name]['part_list']
		for p in part_list:
			cur_part_type = self.parts[p['part_name']]['type']
			if cur_part_type in mapping:
				rev = ''
				if p['dir'] == 'R':
					rev = '<'
				part_col = 13
				if col_mapping != None and p['part_name'] in col_mapping.keys():
					part_col = col_mapping[p['part_name']]
				f_out.write(rev + mapping[cur_part_type] + ' ' + p['part_name'] + ' ' + str(part_col) + '\n')
		f_out.write('# Arcs\n')
		f_out.close()

	def transcriptional_units (self, non_terminated=False):
		"""Extract all valid transcriptional units from a GeneClusterLibrary.

		Note: double promoters will generate two transcriptional units.

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
		p_insts = self.find_part_type_instances('Promoter')
		# Find all the downstream terminators from the promoters
		t_insts = self.find_next_part_idxs(p_insts, part_type='Terminator', dir_check=True)
		# Combine the data to give the full transcriptional unit boundaries
		units = {}
		for v_key in p_insts.keys():
			if v_key in t_insts.keys():
				cur_p_data = p_insts[v_key]
				cur_t_data = t_insts[v_key]
				for idx in range(len(cur_p_data)):
					if cur_t_data[idx] != None or non_terminated == True:
						if v_key not in units.keys():
							units[v_key] = []
						units[v_key].append([cur_p_data[idx], cur_t_data[idx]])
		return units

	def monocistronic_units (self):
		"""Extract all monocistronic transcriptional units from a GeneClusterLibrary.

	    Returns
	    -------
	    units: dict(list([PromoterID, TerminatorID]))
	        All monocistronic transcriptional units. Dictionary indexed by variant name 
	        with elements stored as a list of promoter, terminator pairs (also stored as 
	        a list).
		"""
		return self.polycistronic_units(number_of_CDS=1)

	def polycistronic_units (self, number_of_CDS=None):
		"""Extract all polycistronic transcriptional units from a GeneClusterLibrary.

		Parameters
	    ----------
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
		p_insts = self.find_part_type_instances('Promoter')
		# Find all the downstream terminators from the promoters
		t_insts = self.find_next_part_idxs(p_insts, part_type='Terminator', dir_check=True)
		# For each valid pair, count number of CDS part types between only include if == 1
		units = {}
		for v_key in p_insts.keys():
			if v_key in t_insts.keys():
				cur_p_data = p_insts[v_key]
				cur_t_data = t_insts[v_key]
				for idx in range(len(cur_p_data)):
					cur_p = cur_p_data[idx]
					cur_t = cur_t_data[idx]
					if cur_t != None:
						# Make sure we cycle in the right direction
						step_dir = 1   
						if cur_t < cur_p:
						    step_dir = -1
						# Cycle through all parts between promoter and terminator and count CDSs
						cds_count = 0
						for cur_part_idx in range(cur_p+step_dir, cur_t, step_dir):
							cur_part_name = self.variants[v_key]['part_list'][cur_part_idx]['part_name']
							if self.parts[cur_part_name]['type'] == 'CDS':
								cds_count += 1
						if number_of_CDS == None:
							if cds_count > 1:
								if v_key not in units.keys():
									units[v_key] = []
								units[v_key].append([cur_p, cur_t])
						elif number_of_CDS == cds_count:
							if v_key not in units.keys():
								units[v_key] = []
							units[v_key].append([cur_p, cur_t])
		return units

	def extract_ranges_for_transcriptional_units (self, units):
		"""Extract sequences ranges for a set of transcriptional units.

		Parameters
	    ----------
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
				p_start_idx = self.variants[v_key]['part_list'][el[0]]['seq_idx']
				p_end_idx = p_start_idx + self.variants[v_key]['part_list'][el[0]]['seq_len']
				t_start_idx = self.variants[v_key]['part_list'][el[1]]['seq_idx']
				t_end_idx = t_start_idx + self.variants[v_key]['part_list'][el[1]]['seq_len']
				if v_key not in ranges.keys():
					ranges[v_key] = []
				# We return the full length of the transcriptional unit promoter start -> terminator end
				ranges[v_key].append([p_start_idx, t_end_idx])
		return ranges

	def divergent_promoters (self):
		"""Extract all divergent promoters from a GeneClusterLibrary.

	    Returns
	    -------
	    locations: dict(list(list))
	        Locations of divergent promoters in the library. Stored as dictionary keyed
	        on variant and then a list of tuples containing the locations of the 
	        promoter pairs.
		"""
		results = {}
		# Find all promoters in the library
		p_insts = self.find_part_type_instances('Promoter')
		# For each promoter search upstream to see if another promoter exists
		for v_key in p_insts.keys():
		    for cur_p in p_insts[v_key]:
		        prev_p = self.find_prev_part_idx(v_key, cur_p, part_type='Promoter', next_count=1)
		        if prev_p != None:
		        	# Compare directions
		            cur_p_dir = self.variants[v_key]['part_list'][cur_p]['dir']
		            prev_p_dir = self.variants[v_key]['part_list'][prev_p]['dir']
		            if cur_p_dir != prev_p_dir:
			            step_dir = 1
			            if prev_p < cur_p:
			                step_dir = -1
			            # Check that no terminators, CDSs or RBSs between (ignore spacers and scars)
			            valid_promoter = True
			            for to_check_idx in range(cur_p+step_dir, prev_p, step_dir):
			                part_name = self.variants[v_key]['part_list'][to_check_idx]['part_name']
			                if self.parts[part_name]['type'] in ['CDS', 'RBS', 'Terminator']:
			                    valid_promoter = False
			                    break
			            if valid_promoter == True:
			                # Check to see if already exists (if not then add)
			                if v_key not in results.keys():
			                    results[v_key] = []
			                exists = False
			                for el in results[v_key]:
			                    if el == [cur_p, prev_p] or el == [prev_p, cur_p]:
			                        exists = True
			                        break
			                if exists == False:
			                    results[v_key].append([cur_p, prev_p])
		return results

	def convergent_promoters (self):
		"""Extract all convergent promoters from a GeneClusterLibrary.

	    Returns
	    -------
	    locations: dict(list(list))
	        Locations of convergent promoters in the library. Stored as dictionary keyed
	        on variant and then a list of tuples containing the locations of the 
	        promoter pairs.
		"""
		results = {}
		# Find all promoters in the library
		p_insts = self.find_part_type_instances('Promoter')
		# For each promoter search upstream to see if another promoter exists
		for v_key in p_insts.keys():
		    for cur_p in p_insts[v_key]:
		        next_p = self.find_next_part_idx(v_key, cur_p, part_type='Promoter', next_count=1)
		        if next_p != None:
		            # Compare directions
		            cur_p_dir = self.variants[v_key]['part_list'][cur_p]['dir']
		            next_p_dir = self.variants[v_key]['part_list'][next_p]['dir']
		            if cur_p_dir != next_p_dir:
		                # Promoters are convergent so add.
		                # Check to see if already exists (if not then add)
		                if v_key not in results.keys():
		                    results[v_key] = []
		                exists = False
		                for el in results[v_key]:
		                    if el == [cur_p, next_p] or el == [next_p, cur_p]:
		                        exists = True
		                        break
		                if exists == False:
		                    results[v_key].append([cur_p, next_p])
		return results

	def double_parts (self, part_type):
	    """Extract all double promoters from a GeneClusterLibrary.

	    Returns
	    -------
	    locations: dict(list(list))
	        Locations of double promoters in the library. Stored as dictionary keyed
	        on variant and then a list of tuples containing the locations of the 
	        promoters.
	    """
	    results = {}
	    # Find all promoters in the library
	    p_insts = self.find_part_type_instances(part_type)
	    # For each promoter search upstream to see if another promoter exists
	    for v_key in p_insts.keys():
	        for cur_p in p_insts[v_key]:
	            next_p = self.find_next_part_idx(v_key, cur_p, part_type=part_type, next_count=1)
	            if next_p != None:
	                # Compare directions
	                cur_p_dir = self.variants[v_key]['part_list'][cur_p]['dir']
	                next_p_dir = self.variants[v_key]['part_list'][next_p]['dir']
	                if cur_p_dir == next_p_dir:
	                    # Check that no terminators, CDSs or RBSs between (ignore spacers and scars)
	                    step_dir = 1
	                    if next_p < cur_p:
	                        step_dir = -1
	                    valid_part = True
	                    for to_check_idx in range(cur_p+step_dir, next_p, step_dir):
	                        part_name = self.variants[v_key]['part_list'][to_check_idx]['part_name']
	                        if self.parts[part_name]['type'] in [x for x in ['Promoter', 'CDS', 'RBS', 'Terminator'] if x != part_type]:
	                            valid_part = False
	                            break
	                    if valid_part == True:
	                        # Check to see if already exists (if not then add)
	                        if v_key not in results.keys():
	                            results[v_key] = []
	                        exists = False
	                        for el in results[v_key]:
	                            if el == [cur_p, next_p] or el == [next_p, cur_p]:
	                                exists = True
	                                break
	                        if exists == False:
	                            results[v_key].append([cur_p, next_p])
	    return results

	def double_promoters (self):
	    """Extract all double promoters from a GeneClusterLibrary.

	    Returns
	    -------
	    locations: dict(list(list))
	        Locations of double promoters in the library. Stored as dictionary keyed
	        on variant and then a list of tuples containing the locations of the 
	        promoters.
	    """
	    return self.double_parts('Promoter')

	def double_terminators (self):
		"""Extract all double terminators from a GeneClusterLibrary.

	    Returns
	    -------
	    locations: dict(list(list))
	        Locations of double terminators in the library. Stored as dictionary keyed
	        on variant and then a list of tuples containing the locations of the 
	        terminators.
		"""
		return self.double_parts('Terminator')

	def convergent_terminators (self):
	    """Extract all convergent terminators from a GeneClusterLibrary.

	    Returns
	    -------
	    locations: dict(list(list))
	        Locations of convergent terminators in the library. Stored as dictionary keyed
	        on variant and then a list of tuples containing the locations of the 
	        terminators.
	    """
	    results = {}
	    # Find all promoters in the library
	    p_insts = self.find_part_type_instances('Terminator')
	    # For each promoter search upstream to see if another promoter exists
	    for v_key in p_insts.keys():
	        for cur_p in p_insts[v_key]:
	            next_p = self.find_next_part_idx(v_key, cur_p, part_type='Terminator', next_count=1)
	            if next_p != None:
	                # Compare directions
	                cur_p_dir = self.variants[v_key]['part_list'][cur_p]['dir']
	                next_p_dir = self.variants[v_key]['part_list'][next_p]['dir']
	                if cur_p_dir != next_p_dir:
	                    # Terminators are convergent so add.
	                    # Check that no terminators, CDSs or RBSs between (ignore spacers and scars)
	                    step_dir = 1
	                    if next_p < cur_p:
	                        step_dir = -1
	                    valid_part = True
	                    for to_check_idx in range(cur_p+step_dir, next_p, step_dir):
	                        part_name = self.variants[v_key]['part_list'][to_check_idx]['part_name']
	                        if self.parts[part_name]['type'] in ['Promoter', 'CDS', 'RBS']:
	                            valid_part = False
	                            break
	                    if valid_part == True:
		                    # Check to see if already exists (if not then add)
		                    if v_key not in results.keys():
		                        results[v_key] = []
		                    exists = False
		                    for el in results[v_key]:
		                        if el == [cur_p, next_p] or el == [next_p, cur_p]:
		                            exists = True
		                            break
		                    if exists == False:
		                        results[v_key].append([cur_p, next_p])
	    return results

	def tu_meta_data (self, variant, tu):
		"""Find meta data for the transcriptional unit.

		If any of the following components are not found then None is returned.

		Parameters
		----------
		variant : string
			Gene cluster variant name.

		tu : list([start_part_idx, end_part_idx])
			Transcriptional unit start and end part indexes.

		Returns
		-------
		metadata : list(string) ([promoter, rbs, cds, terminator])
			Metadata for the transcriptional unit (part names)
		"""
		p = None
		r = None
		c = None
		t = None
		# Make sure we cycle in the right direction
		step_dir = 1
		if tu[1] < tu[0]:
		    step_dir = -1
		# Cycle through all parts to extract meta data
		for cur_part_idx in range(tu[0], tu[1]+step_dir, step_dir):
			cur_name = self.variants[variant]['part_list'][cur_part_idx]['part_name']
			cur_type = self.parts[cur_name]['type']
			# Check to see if found relevant part
			if p == None and cur_type == 'Promoter':
				p = cur_name
			if r == None and cur_type == 'RBS':
				r = cur_name
			if c == None and cur_type == 'CDS':
				c = cur_name
			if t == None and cur_type == 'Terminator':
				t = cur_name
		return [p, r, c, t]

	def tu_insts_meta_data (self, tu_insts):
		"""Find meta data for a set of transcriptional unit instances.

		If any of the following components are not found then None is returned.

		Parameters
		----------
		variant : string
			Gene cluster variant name.

		tu : list([start_part_idx, end_part_idx])
			Transcriptional unit start and end part indexes.

		Returns
		-------
		metadata : dict(list(string)) ([promoter, rbs, cds, terminator])
			Metadata for the transcriptional units (part names) with dictionary
			keyed by variant name.
		"""
		results = {}
		for v_key in tu_insts.keys():
			results[v_key] = []
			for tu in tu_insts[v_key]:
				results[v_key].append(self.tu_meta_data(v_key, tu))
		return results

	def __make_float_if_needed (self, s):
		"""Helper function to automatically convert a string to a number if possible.

	    Parameters
	    ----------
		s : string
			Value to attempt conversion to float.

	    Returns
	    -------
	    s : string or float
	    	If successful the float version will be returned, otherwise the original
	    	string is returned.
		"""
		try:
			float(s)
			return float(s)
		except ValueError:
			return s

	def __remove_nones_variant_data (self, variant_data):
		"""Helper function to remove None values from variant data.

		Many of the internal functions return results as a dictionary separating
		them into the specific variants in which they are found. Sometimes these
		will include None value elements. This function will remove these elements.

	    Parameters
	    ----------
		variant_data : dict(list)
			Variant data to remove None values from.

	    Returns
	    -------
	    variant_data : dict(list)
	    	Cleaned version of variant_data with no None values.
		"""
		for v in variant_data.keys():
			variant_data[v] = [x for x in variant_data[v] if x is not None]
		return variant_data

	def __collate_variant_data_list (self, variant_data):
		"""Helper function to collate variant separated data.

		Many of the internal functions return results as a dictionary separating
		them into the specific variants in which they are found. This function 
		will take a dictionary of data values and collate all the data lists 
		into a single list.

	    Parameters
	    ----------
		variant_data : dict(list)
			Variant data to collate.

	    Returns
	    -------
	    variant_data : list
	    	Single list containing the collated data.
		"""
		collated = []
		for v in variant_data.keys():
			collated = collated + variant_data[v]
		return collated

	def __collate_variant_data_value (self, variant_data):
		"""Helper function to collate variant separated data.

		Many of the internal functions return results as a dictionary separating
		them into the specific variants in which they are found. This function 
		will take a dictionary of data values and collate all the data values 
		into a single list.

	    Parameters
	    ----------
		variant_data : dict(list)
			Variant data to collate.

	    Returns
	    -------
	    variant_data : list
	    	Single list containing the collated data.
		"""
		collated = []
		for v in variant_data.keys():
			collated.append(variant_data[v])
		return collated

###############################################################################
# GLOBAL FUNCTIONS
###############################################################################

def __make_float_if_needed (s):
	"""Helper function to automatically convert a string to a number if possible.

    Parameters
    ----------
	s : string
		Value to attempt conversion to float.

    Returns
    -------
    s : string or float
    	If successful the float version will be returned, otherwise the original
    	string is returned.
	"""
	try:
		float(s)
		return float(s)
	except ValueError:
		return s

def load_part_attribs (parts_filename):
	"""Load part attributes from CSV file.

    Parameters
    ----------
	parts_filename : string
		CSV filename.

    Returns
    -------
    part_data : dict(list)
    	Converts the CSV file into a dictionary of parts with each key returning a
    	list of the part details and attributes.
	"""
	part_data = {}
	parts_reader = csv.reader(open(parts_filename, 'rU'), delimiter=',')
	header = next(parts_reader)
	header_map = {}
	for i in range(len(header)):
		header_map[header[i]] = i
	attrib_keys = [k for k in header_map.keys() if k not in ['Part', 'Part Type', 'Part Sequence']]
	for row in parts_reader:
		# Make the attributes map
		part_attribs_map = {}
		for k in attrib_keys:
			if row[header_map[k]] != '':
				part_attribs_map[k] = __make_float_if_needed(row[header_map[k]])
		part_name = row[header_map['Part']]
		part_type = row[header_map['Part Type']]
		part_seq = row[header_map['Part Sequence']]
		part_data[part_name] = [part_type, part_seq, part_attribs_map]
	return part_data

def load_variant_attribs (variant_attribs_filename):
	"""Load variant attributes from CSV file.

    Parameters
    ----------
	variant_attribs_filename : string
		CSV filename.

    Returns
    -------
    variant_attribs_data : dict(dict)
    	Converts the CSV file into a dictionary of variants with each key returning a
    	list of the variant attributes also stored as a dictionary.
	"""
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
			if row[header_map[k]] != '':
				variant_attribs_map[k] = __make_float_if_needed(row[header_map[k]])
		variant_name = row[header_map['Variant']]
		variant_attribs_data[variant_name] = variant_attribs_map
	return variant_attribs_data

def load_variant_designs (variant_designs_filename):
	"""Load variant designs from CSV file.

    Parameters
    ----------
	variant_designs_filename : string
		CSV filename.

    Returns
    -------
    variant_designs_data : dict(list)
    	Converts the CSV file into a dictionary of variants with each key returning a
    	list of the part names that make up the variant design.
	"""
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
	"""Load variant part attributes from CSV file.

    Parameters
    ----------
	variant_part_attribs_filename : string
		CSV filename.

    Returns
    -------
    variant_part_attribs_data : dict(dict(dict(string)))
    	Converts the CSV file into a dictionary of dictionaries with the keys:
    	variant, part, attribute and the element returned the attribute value. This
    	allows for data to be stored at the variant-part level e.g., transcriptomic
    	data for genes.
	"""
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
			if row[header_map[part]] != '':
				variant_part_attribs_data[variant_name][part][attrib_name] = __make_float_if_needed(row[header_map[part]])
	return variant_part_attribs_data

def make_library_from_csv (part_attribs_fn, variant_designs_fn, variant_attribs_fn, 
	                       variant_part_attribs_fn, output_fn):
	"""Create a GeneClusterLibrary from a set of CSV files.

    Parameters
    ----------
	part_attribs_fn : string
		Part attributes CSV filename.

	variant_designs_fn : string
		Variant designs CSV filename.

	variant_attribs_fn : string
		Variant attributes CSV filename.

	variant_part_attribs_fn : string
		Variant part attributes CSV filename.

	output_fn : string
		Output file containing the GeneClusterLibrary.
	"""
	# Load the data
	part_attribs = load_part_attribs(part_attribs_fn)
	variant_designs = load_variant_designs(variant_designs_fn)
	variant_attribs = load_variant_attribs(variant_attribs_fn)
	variant_part_attribs = load_variant_part_attribs(variant_part_attribs_fn)
	# Create GeneClusterLibrary object
	gene_cluster_lib = GeneClusterLibrary()
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
