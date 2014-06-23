#!/usr/bin/env python
"""
Gene Cluster Tools
==================

    Gene Cluster Tools is a small library designed for creating, querying
    and saving information related to refactored gene clusters. It provides
    data structures and functions to manipulate light-weight text based
    files that contain all sequence, part and cluster variant information,
    in addition to search capabilities to facilitate easy querying of data
    during analysis. It can output clusters to PigeonCAD format to generate
    a visual representation.
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

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

class GeneClusterLibrary:
	"""Class encapsulating the concept of a gene cluster and potential variants.
	"""

	# Static class variable holding standard mapping of part type to PigeonCAD type
	STD_PIGEON_MAP = {'Promoter':'p',
	                  'CDS':'c',
	                  'RBS':'r',
	                  'Terminator':'t',
	                  'Scar':'x',
	                  'Spacer':'s'}

	def __init__(self):
		"""Constructor the generate an empty GeneClusterLibrary.
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
		of a name, type and sequence, and can additional attributes associated
		with it e.g., performance information such as REU. All attributes are
		held as a dict to enable easy and fast keyed access. Parts that have
		been defined can be accessed using the 'parts' instance variable.

	    Parameters
	    ----------
	    name : string
	        Part name.

	    part_type : string
	    	Part type (see STD_PIGEON_MAP for examples)

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
	    seqs : dict(list(string))
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

	def find_next_part_idx (self, variant_name, part_idx, part_type=None, next_count=1):
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
						if cur_next_count == next_count:
							return i
						else:
							cur_next_count += 1
		return None

	def find_next_part_idxs (self, variant_insts, part_type=None, next_count=1, remove_nones=False):
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
				v_idxs.append(self.find_next_part_idx(v, i, part_type, next_count))
			part_idxs[v] = v_idxs
		if remove_nones:
			return self.__remove_nones_variant_data(part_idxs)
		else:
			return part_idxs

	def find_prev_part_idx (self, variant_name, part_idx, part_type=None, next_count=1):
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
						if cur_next_count == next_count:
							return i
						else:
							cur_next_count += 1
		return None
	
	def find_prev_part_idxs (self, variant_insts, part_type=None, next_count=1, remove_nones=False):
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
				v_idxs.append(self.find_prev_part_idx(v, i, part_type, next_count))
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

	def find_next_parts (self, variant_insts, n_next=1):
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
			part_insts = self.find_next_part_idxs(variant_insts, next_count=i, remove_nones=True)
			next_parts = next_parts + self.extract_part_names_from_idxs(part_insts)
		return next_parts

	def find_prev_parts (self, variant_insts, n_prev=1):
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
			part_insts = self.find_prev_part_idxs(variant_insts, next_count=i, remove_nones=True)
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

	def __make_float_if_needed(self, s):
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



def promoter_strength (trace, edge_start=50, edge_skip_end=200, edge_end=50):
	"""Calculate the promoter strength (delta F) of trace
	"""
	start_depth = np.median(trace[0:edge_start])
	end_depth = np.median(trace[-edge_length:])
	return None



