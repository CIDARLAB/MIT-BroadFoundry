#!/usr/bin/env python
"""
Analysis of the raw RNA-seq data from the circuits group.

Conditions:
	t1 = none
	t2 = iptg                 (PTac)
	t3 = atc                  (PTet)
	t4 = iptg + atc           (PTac, PTet)
	t5 = ara                  (PBAD)
	t6 = ara + iptg           (PBAD, PTac)
	t7 = ara + atc            (PBAD, PTet)
	t8 = ara + iptg + atc     (PBAD, PTac, PTet)
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv

matplotlib.rcParams['lines.dash_joinstyle']  = 'miter'
matplotlib.rcParams['lines.dash_capstyle']   = 'butt'
matplotlib.rcParams['lines.solid_joinstyle'] = 'miter'
matplotlib.rcParams['lines.solid_capstyle']  = 'projecting'
# Make text editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype']          = 42

DATA_PREFIX = '../results'
OUTPUT_PREFIX = './plots_transfer'

MIN_READS_TERM = 500.0
MIN_READS_RIBO = 500.0

tube_samples = ['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8']
flask_samples = ['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8']

plt.rcParams['ytick.major.pad']='1'
plt.rcParams['xtick.major.pad']='5'
fmt_label_size = 13.5
fmt_edge_width = 2.5
trace_height = 2.5

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

def load_fpkm_matrix (filename):
	fpkm_data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	fpkm_idx = {}
	for idx in range(len(header)-1):
		fpkm_idx[header[idx+1]] = idx+1
	for row in data_reader:
		if len(row) == 17:
			fpkm_data[row[0]] = [float(x) for x in row[1:]]
	return fpkm_data, fpkm_idx

def plot_terminator_transfer (data, samples, fpkm_in, filename_out, point_shape='o', min_upstream_read=None):
	plt.rcParams['ytick.major.pad']='1'
	plt.rcParams['xtick.major.pad']='5'
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	y_vals = []
	for data_point in data:
		if data_point[0] in samples:
			d_us_reads = float(data_point[1])
			d_ds_reads = float(data_point[2])
			d_t_e = float(data_point[3])
			d_t_s = float(data_point[4])
			d_max_term = data_point[5]
			good_val = True
			point_col = (0,0,0)
			if d_us_reads < min_upstream_read:
				point_col = (0.7,0.7,0.7)
				good_val = False
			elif d_max_term == 'Y':
				point_col = (0.95,0.30,0.25)
			if good_val == True:
				y_vals.append(d_t_s)
			ax.scatter([fpkm_in[data_point[0]]], [d_t_s], marker=point_shape, zorder=100, s=40, color=(1,1,1,0), edgecolor=point_col, linewidth=2.0)
	plt.axhline(y=np.median(y_vals), linestyle='--', color=(0,0,0), linewidth=fmt_char_data_linewidth)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([30,200000])
	ax.set_xscale('log')
	ax.set_ylim([0.5,20000])
	ax.set_yscale('symlog', linthreshy=1.0)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_edge_width)
		ax.spines[axis].set_linewidth(fmt_edge_width)
		ax.spines[axis].set_linewidth(fmt_edge_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	#ax.set_xticks([])
	#ax.set_xticklabels([], visible=False)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')

def plot_ribozyme_transfer (data, samples, fpkm_out, filename_out, point_shape='o', min_downstream_read=None):
	plt.rcParams['ytick.major.pad']='5'
	plt.rcParams['xtick.major.pad']='5'
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	x_vals = []
	y_vals = []
	for data_point in data:
		if data_point[0] in samples:
			d_us_reads = float(data_point[1])
			d_ds_reads = float(data_point[2])
			d_c_e = float(data_point[3])
			d_max_term = data_point[4]
			d_cut_site = int(data_point[5])
			good_val = True
			point_col = (0,0,0)
			if d_ds_reads < min_downstream_read:
				point_col = (0.7,0.7,0.7)
				good_val = False
			elif d_max_term == 'Y':
				point_col = (0.95,0.30,0.25)
			if good_val == True:
				x_vals.append(d_ds_reads)
				y_vals.append(d_c_e)
			ax.scatter([fpkm_out[data_point[0]]], [d_c_e], marker=point_shape, zorder=100, s=40, color=(1,1,1,0), edgecolor=point_col, linewidth=2.0)
	plt.axhline(y=np.median(y_vals), linestyle='--', color=(0,0,0), linewidth=fmt_char_data_linewidth)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([30,200000])
	ax.set_xscale('symlog')
	ax.set_ylim([-0.1,1.1])
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_edge_width)
		ax.spines[axis].set_linewidth(fmt_edge_width)
		ax.spines[axis].set_linewidth(fmt_edge_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')


def list_to_dict(the_list, dict_keys):
	new_dict = {}
	for el_idx in range(len(the_list)):
		new_dict[dict_keys[el_idx]] = the_list[el_idx]
	return new_dict

fpkm_data, fpkm_idx = load_fpkm_matrix(DATA_PREFIX+'/fpkm.normed.matrix.txt')

########### TERMINATORS ###########
terminator_data = load_data(DATA_PREFIX+'/terminator.profile.perf.txt')
terminators = ['ECK120029600', 'ECK120033737', 'L3S2P11', 'L3S2P21', 'L3S2P24', 'L3S2P55', 'L3S2P22', 'L3S3P21-2']
term_tube_fpkms_in = {'ECK120029600': list_to_dict(fpkm_data['SrpR'][8:16],tube_samples),
                      'ECK120033737': list_to_dict(fpkm_data['PhlF'][8:16],tube_samples),
                      'L3S2P11': list_to_dict(fpkm_data['BM3R1'][8:16],tube_samples),
                      'L3S2P21': list_to_dict(fpkm_data['YFP'][8:16],tube_samples),
                      'L3S2P24': list_to_dict(fpkm_data['LitR'][8:16],tube_samples),
                      'L3S2P55': list_to_dict(fpkm_data['AmtR'][8:16],tube_samples),
                      'L3S2P22': list_to_dict(fpkm_data['AraC'][8:16],tube_samples),
                      'L3S3P21-2': list_to_dict(fpkm_data['KanR'][8:16],tube_samples)}
term_flask_fpkms_in = {'ECK120029600': list_to_dict(fpkm_data['SrpR'][0:8],flask_samples),
                      'ECK120033737': list_to_dict(fpkm_data['PhlF'][0:8],flask_samples),
                      'L3S2P11': list_to_dict(fpkm_data['BM3R1'][0:8],flask_samples),
                      'L3S2P21': list_to_dict(fpkm_data['YFP'][0:8],flask_samples),
                      'L3S2P24': list_to_dict(fpkm_data['LitR'][0:8],flask_samples),
                      'L3S2P55': list_to_dict(fpkm_data['AmtR'][0:8],flask_samples),
                      'L3S2P22': list_to_dict(fpkm_data['AraC'][0:8],flask_samples),
                      'L3S3P21-2': list_to_dict(fpkm_data['KanR'][0:8],flask_samples)}
for t in terminators:
	plot_terminator_transfer(terminator_data['0x58v50'][t], tube_samples, term_tube_fpkms_in[t], OUTPUT_PREFIX+'/term_'+t+'_tube.pdf', point_shape='o', min_upstream_read=MIN_READS_TERM)
	plot_terminator_transfer(terminator_data['0x58v50'][t], flask_samples, term_flask_fpkms_in[t], OUTPUT_PREFIX+'/term_'+t+'_flask.pdf', point_shape='x', min_upstream_read=MIN_READS_TERM)

########### RIBOZYMES ###########
ribozyme_data = load_data(DATA_PREFIX+'/ribozyme.profile.perf.txt')
ribozymes = ['BydvJ', 'PlmJ', 'SarJ', 'RiboJ10', 'RiboJ53', 'RiboJ']
ribo_tube_fpkms_out = {'BydvJ': list_to_dict(fpkm_data['AmtR'][8:16],tube_samples),
                      'PlmJ': list_to_dict(fpkm_data['LitR'][8:16],tube_samples),
                      'SarJ': list_to_dict(fpkm_data['BM3R1'][8:16],tube_samples),
                      'RiboJ10': list_to_dict(fpkm_data['SrpR'][8:16],tube_samples),
                      'RiboJ53': list_to_dict(fpkm_data['PhlF'][8:16],tube_samples),
                      'RiboJ': list_to_dict(fpkm_data['YFP'][8:16],tube_samples)}
ribo_flask_fpkms_out = {'BydvJ': list_to_dict(fpkm_data['AmtR'][0:8],flask_samples),
                      'PlmJ': list_to_dict(fpkm_data['LitR'][0:8],flask_samples),
                      'SarJ': list_to_dict(fpkm_data['BM3R1'][0:8],flask_samples),
                      'RiboJ10': list_to_dict(fpkm_data['SrpR'][0:8],flask_samples),
                      'RiboJ53': list_to_dict(fpkm_data['PhlF'][0:8],flask_samples),
                      'RiboJ': list_to_dict(fpkm_data['YFP'][0:8],flask_samples)}
for r in ribozymes:
	plot_ribozyme_transfer(ribozyme_data['0x58v50'][r], tube_samples, ribo_tube_fpkms_out[r], OUTPUT_PREFIX+'/ribo_'+r+'_tube.pdf', point_shape='o', min_downstream_read=MIN_READS_RIBO)
	plot_ribozyme_transfer(ribozyme_data['0x58v50'][r], flask_samples, ribo_flask_fpkms_out[r], OUTPUT_PREFIX+'/ribo_'+r+'_flask.pdf', point_shape='x', min_downstream_read=MIN_READS_RIBO)


########### ROBUSTNESS DATA ###########

gene_to_term = {}
gene_to_term['AmtR'] = 'L3S2P55'
gene_to_term['LitR'] = 'L3S2P24'
gene_to_term['BM3R1'] = 'L3S2P11'
gene_to_term['SrpR'] = 'ECK120029600'
gene_to_term['PhlF'] = 'ECK120033737'
gene_to_term['YFP'] = 'L3S2P21'

f_out_1 = open(OUTPUT_PREFIX+'/term_data_1.txt', 'w')
f_out_2 = open(OUTPUT_PREFIX+'/term_data_2.txt', 'w')
for g in gene_to_term.keys():
	t_e_list = terminator_data['0x58v50'][gene_to_term[g]]
	for el in t_e_list:
		t_e = str(1.0-float(el[3]))
		if el[0] in tube_samples:
			f_out_1.write(' '.join([g,'t',t_e]) + '\n')
		else:
			f_out_2.write(' '.join([g,'t',t_e]) + '\n')
f_out_1.close()
f_out_2.close()



