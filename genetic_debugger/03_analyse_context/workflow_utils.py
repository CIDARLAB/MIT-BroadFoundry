#!/usr/bin/env python
"""Utils for assessing contextual effects
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def filter_contexts (p_data, context_list):
	new_data = []
	for i in p_data:
		if i[3] in context_list:
			new_data.append(i)
	return new_data

def filter_part_names (p_data, part_name_list):
	new_data = []
	for i in p_data:
		if i[2] in part_name_list:
			new_data.append(i)
	return new_data

def filter_details (p_data, key, detail_list):
	new_data = []
	for i in p_data:
		if i[4][key] in detail_list:
			new_data.append(i)
	return new_data

def all_part_names (p_data):
	return sorted(list(set([x[2] for x in p_data])))

def all_contexts (p_data):
	return sorted(list(set([x[3] for x in p_data])))

def all_details (p_data, key):
	return sorted(list(set([x[4][key] for x in p_data if x[4][key] != None])))

def split_on_context (p_data):
	context_data = {}
	for context in all_contexts(p_data):
		context_data[context] = filter_contexts(p_data, [context])
	return context_data

def split_on_part (p_data):
	part_data = {}
	for part in all_part_names(p_data):
		part_data[part] = filter_part_names(p_data, [part])
	return part_data

def split_on_detail (p_data, key):
	detail_data = {}
	for detail in all_details(p_data, key):
		detail_data[detail] = filter_details(p_data, key, [detail])
	return detail_data

def all_in_fluxes (p_data):
	return [x[5][0] for x in p_data]

def all_out_fluxes (p_data):
	return [x[5][1] for x in p_data]

def all_strengths (p_data):
	return [x[5][2] for x in p_data]

def assess_contextual_effects (p_data):

	# Overall variability in expression across all promoter types
	
	# Strengths break down by promoter type
	
	# How are these effected by architecture
	
	# How are these effected by up and downstream elements
	
	# Do architecture and elements up/downstream work synegistically
	
	return None




