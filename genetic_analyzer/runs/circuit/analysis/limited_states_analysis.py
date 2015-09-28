#!/usr/bin/env python
"""
Analysis of part performance from limited numbers of input states.
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
import math
import csv

DATA_PREFIX = '../results'
OUTPUT_PREFIX = './limited_states'

MIN_READS_TERM = 500.0
MIN_READS_RIBO = 500.0

tube_samples = ['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8']
flask_samples = ['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8']

def load_data (filename):
	data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) > 3:
			chrom = row[0]
			part_name = row[1]
			if chrom not in data.keys():
				data[chrom] = {}
			if part_name not in data[chrom].keys():
				data[chrom][part_name] = []
			data[chrom][part_name].append(row[2:])
	return data

def calc_terminator_perf (data, samples, min_upstream_read=None):
	y_vals = []
	for data_point in data:
		if data_point[0] in samples:
			d_us_reads = float(data_point[1])
			d_ds_reads = float(data_point[2])
			d_t_s = float(data_point[4])
			if d_us_reads > min_upstream_read:
				y_vals.append(d_t_s)
	if y_vals != []:
		return np.median(y_vals)
	else:
		return None

def calc_ribozyme_perf (data, samples, min_downstream_read=None):
	y_vals = []
	for data_point in data:
		if data_point[0] in samples:
			d_us_reads = float(data_point[1])
			d_ds_reads = float(data_point[2])
			d_c_e = float(data_point[3])
			if d_ds_reads > min_downstream_read:
				y_vals.append(d_c_e)
	if y_vals != []:
		return np.median(y_vals)
	else:
		return None

# Sets of states to calculate performance for
sub_states_tube  = [['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8'],
                    ['tube_1', 'tube_2', 'tube_3', 'tube_5', 'tube_6'],
                    ['tube_1', 'tube_2', 'tube_4', 'tube_5', 'tube_6'],
                    ['tube_1', 'tube_2', 'tube_3', 'tube_5', 'tube_7'],
                    ['tube_1', 'tube_2', 'tube_4', 'tube_5', 'tube_7'],
                    ['tube_1', 'tube_2', 'tube_3', 'tube_5', 'tube_8'],
                    ['tube_1', 'tube_2', 'tube_4', 'tube_5', 'tube_8']]

sub_states_flask = [['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8'],
                    ['flask_1', 'flask_2', 'flask_3', 'flask_5', 'flask_6'],
                    ['flask_1', 'flask_2', 'flask_4', 'flask_5', 'flask_6'],
                    ['flask_1', 'flask_2', 'flask_3', 'flask_5', 'flask_7'],
                    ['flask_1', 'flask_2', 'flask_4', 'flask_5', 'flask_7'],
                    ['flask_1', 'flask_2', 'flask_3', 'flask_5', 'flask_8'],
                    ['flask_1', 'flask_2', 'flask_4', 'flask_5', 'flask_8']]

########### TERMINATORS ###########
terminator_data = load_data(DATA_PREFIX+'/terminator.profile.perf.txt')
terminators = ['ECK120029600', 'ECK120033737', 'L3S2P11', 'L3S2P21', 'L3S2P24', 'L3S2P55', 'L3S2P22', 'L3S3P21-2']
terminator_tube_ss_perf = {}
terminator_flask_ss_perf = {}
for t in terminators:
	term_perfs_tube = []
	for ss in sub_states_tube:
		term_perfs_tube.append( calc_terminator_perf(terminator_data['0x58v50'][t], ss, min_upstream_read=MIN_READS_TERM) )
	terminator_tube_ss_perf[t] = term_perfs_tube
	term_perfs_flask = []
	for ss in sub_states_flask:
		term_perfs_flask.append( calc_terminator_perf(terminator_data['0x58v50'][t], ss, min_upstream_read=MIN_READS_TERM) )
	terminator_flask_ss_perf[t] = term_perfs_flask

########### RIBOZYMES ###########
ribozyme_data = load_data(DATA_PREFIX+'/ribozyme.profile.perf.txt')
ribozymes = ['BydvJ', 'PlmJ', 'SarJ', 'RiboJ10', 'RiboJ53', 'RiboJ']
ribozyme_tube_ss_perf = {}
ribozyme_flask_ss_perf = {}
for r in ribozymes:
	ribo_perfs_tube = []
	for ss in sub_states_tube:
		ribo_perfs_tube.append( calc_ribozyme_perf(ribozyme_data['0x58v50'][r], ss, min_downstream_read=MIN_READS_RIBO) )
	ribozyme_tube_ss_perf[r] = ribo_perfs_tube
	ribo_perfs_flask = []
	for ss in sub_states_flask:
		ribo_perfs_flask.append( calc_ribozyme_perf(ribozyme_data['0x58v50'][r], ss, min_downstream_read=MIN_READS_RIBO) )
	ribozyme_flask_ss_perf[r] = ribo_perfs_flask

def save_perf_to_file (data, sub_states, filename):
	f_out = open(filename, 'w')
	header = ['part_name']+[','.join(x) for x in sub_states]
	f_out.write('\t'.join(header)+'\n')
	for part_name in data.keys():
		part_data = [str(x) for x in data[part_name]]
		f_out.write('\t'.join([part_name]+part_data)+'\n')
	f_out.close()

save_perf_to_file(terminator_tube_ss_perf, sub_states_tube, OUTPUT_PREFIX+'/term_tube_limited_perf.txt')
save_perf_to_file(terminator_flask_ss_perf, sub_states_flask, OUTPUT_PREFIX+'/term_flask_limited_perf.txt')

save_perf_to_file(ribozyme_tube_ss_perf, sub_states_tube, OUTPUT_PREFIX+'/ribo_tube_limited_perf.txt')
save_perf_to_file(ribozyme_flask_ss_perf, sub_states_flask, OUTPUT_PREFIX+'/ribo_flask_limited_perf.txt')

