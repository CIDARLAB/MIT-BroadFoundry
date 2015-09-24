#!/usr/bin/env python
"""
	Plot transcription profile plots for paper.
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

DATA_PREFIX = '../data'
RESULTS_PREFIX = '../results'
OUTPUT_PREFIX = './plots_profile'

MIN_READS_TERM = 500.0
MIN_READS_RIBO = 500.0

tube_samples = ['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8']
flask_samples = ['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8']

plt.rcParams['ytick.major.pad']='1' # 5 for all but repressor graphs
plt.rcParams['xtick.major.pad']='5'
fmt_label_size = 13.5
fmt_edge_width = 2.5
trace_height = 2.5

def load_gff (filename):
	gff = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process each line
	for row in data_reader:
		if len(row) == 9:
			chromo = row[0]
			part_type = row[2]
			start_bp = int(row[3])
			end_bp = int(row[4])
			part_dir = row[6]
			part_attribs = {}
			split_attribs = row[8].split(';')
			part_name = None
			for attrib in split_attribs:
				key_value = attrib.split('=')
				if len(key_value) == 2:
					if key_value[0] == 'Name':
						part_name = key_value[1]
					else:
						part_attribs[key_value[0]] = key_value[1]
			if part_name != None:
				if chromo not in gff.keys():
					gff[chromo] = {}
				gff[chromo][part_name] = [part_type, part_dir, start_bp, end_bp, part_attribs]
	return gff	

def find_profile (profiles, chrom, start_bp, end_bp):
	if chrom in profiles.keys():
		for el in profiles[chrom]:
			if el[0] == start_bp and el[1] == end_bp:
				return [el[2], el[3]]
	return None

def load_profiles (fwd_profile_filename, rev_profile_filename):
	""" Profiles have the form of a list chr: [start_bp, end_bp, [profile_fwd],[profile_rev]]
	"""
	profiles = {}
	data_reader = csv.reader(open(fwd_profile_filename, 'rU'), delimiter='\t')
	for row in data_reader:
		if len(row) == 5:
			cur_chrom = row[0]
			if cur_chrom not in profiles.keys():
				profiles[cur_chrom] = []
			cur_start_bp = int(row[1])
			cur_end_bp = int(row[2])
			cur_profile = find_profile(profiles, cur_chrom, cur_start_bp, cur_end_bp)
			if cur_profile == None:
				new_profile = [cur_start_bp, cur_end_bp, np.zeros(cur_end_bp-cur_start_bp), np.zeros(cur_end_bp-cur_start_bp)]
				new_profile[2][int(row[3])-1] = int(row[4])
				profiles[cur_chrom].append(new_profile)
			else:
				cur_profile[0][int(row[3])-1] = int(row[4])
	data_reader = csv.reader(open(rev_profile_filename, 'rU'), delimiter='\t')
	for row in data_reader:
		if len(row) == 5:
			cur_chrom = row[0]
			if cur_chrom not in profiles.keys():
				profiles[cur_chrom] = []
			cur_start_bp = int(row[1])
			cur_end_bp = int(row[2])
			cur_profile = find_profile(profiles, cur_chrom, cur_start_bp, cur_end_bp)
			if cur_profile != None:
				cur_profile[1][int(row[3])-1] = int(row[4])
	return profiles

def extract_profile_region (profiles, chrom, start_bp, end_bp):
	region = None
	if chrom in profiles.keys():
		for profile in profiles[chrom]:
			full_chrom = False
			if profile[0] == 0 and profile[1] == len(profile[2]):
				full_chrom = True
			if full_chrom == True:
				fwd_profile = list(profile[2])
				rev_profile = list(profile[3])
				profile_len = len(fwd_profile)
				ext_start_fwd = []
				ext_end_fwd = []
				ext_start_rev = []
				ext_end_rev = []
				# The region will exist
				if start_bp < 0:
					# extend the profile at start
					ext_start_fwd = fwd_profile[start_bp:]
					ext_start_rev = rev_profile[start_bp:]
				if end_bp > profile_len:
					# extend the profile at end
					ext_end_fwd = fwd_profile[:(end_bp-profile_len)]
					ext_end_rev = rev_profile[:(end_bp-profile_len)]
				new_start_bp = start_bp
				new_end_bp = end_bp
				if ext_start_fwd != []:
					new_start_bp = 0
					new_end_bp = end_bp+len(ext_start_fwd)
				new_fwd_profile = ext_start_fwd+fwd_profile+ext_end_fwd
				new_rev_profile = ext_start_rev+rev_profile+ext_end_rev
				region = [new_fwd_profile[new_start_bp:new_end_bp], 
				          new_rev_profile[new_start_bp:new_end_bp]]
				break
	return region

def reverse_region (region):
	return [region[1][::-1], region[0][::-1]]

def plot_terminator_profile (profiles, trace_region, part_region, filename_out):
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	for s in profiles:
		profile_data = extract_profile_region(profiles[s], trace_region[0], trace_region[1], trace_region[2])
		plt.axvspan(part_region[1], part_region[2], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
		ax.plot(range(trace_region[1], trace_region[2]), np.array(profile_data[0]), color=(0.0,0.0,0.0), alpha=1, linewidth=1.5)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([trace_region[1],trace_region[2]])
	ax.set_xscale('linear')
	ax.set_ylim([0,1000000])
	ax.set_yscale('symlog', linthreshy=10)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xticks([])
	ax.set_xticklabels([], visible=False)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')

def plot_ribozyme_profile (profiles, trace_region, part_region, cut_site, filename_out):
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	for s in profiles:
		profile_data = extract_profile_region(profiles[s], trace_region[0], trace_region[1], trace_region[2])
		plt.axvspan(part_region[1], part_region[2], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
		ax.plot(range(trace_region[1], trace_region[2]), np.array(profile_data[0]), color=(0.0,0.0,0.0), alpha=1, linewidth=1.5)
	plt.axvline(x=cut_site, linewidth=1.5, linestyle='--', color=(0,0,0))
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([trace_region[1],trace_region[2]])
	ax.set_xscale('linear')
	ax.set_ylim([0,1000000])
	ax.set_yscale('symlog', linthreshy=10)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xticks([])
	ax.set_xticklabels([], visible=False)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')

def plot_promoter_profile (profiles, trace_region, part_region, filename_out):
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	for s in profiles:
		profile_data = extract_profile_region(profiles[s], trace_region[0], trace_region[1], trace_region[2])
		plt.axvspan(part_region[1], part_region[2], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
		ax.plot(range(trace_region[1], trace_region[2]), np.array(profile_data[0]), color=(0.0,0.0,0.0), alpha=1, linewidth=1.5)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([trace_region[1],trace_region[2]])
	ax.set_xscale('linear')
	ax.set_ylim([0,1000000])
	ax.set_yscale('symlog', linthreshy=10)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xticks([])
	ax.set_xticklabels([], visible=False)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')

# Load design information and transcription profiles
gff = load_gff (DATA_PREFIX+'/gff/0x58v50.gff')
tube_profiles = {}
for s in tube_samples:
	tube_profiles[s] = load_profiles(RESULTS_PREFIX+'/'+s+'/'+s+'.fwd.profiles.txt', RESULTS_PREFIX+'/'+s+'/'+s+'.rev.profiles.txt')
flask_profiles = {}
for s in flask_samples:
	flask_profiles[s] = load_profiles(RESULTS_PREFIX+'/'+s+'/'+s+'.fwd.profiles.txt', RESULTS_PREFIX+'/'+s+'/'+s+'.rev.profiles.txt')

########### TERMINATORS ###########
terminators = ['ECK120029600', 'ECK120033737', 'L3S2P11', 'L3S2P21', 'L3S2P24', 'L3S2P55', 'L3S2P22', 'L3S3P21-2']
for t in terminators:
	trace_region = ['0x58v50', gff['0x58v50'][t][2]-200, gff['0x58v50'][t][3]+50]
	part_region = ['0x58v50', gff['0x58v50'][t][2], gff['0x58v50'][t][3]]
	plot_terminator_profile(tube_profiles, trace_region, part_region, OUTPUT_PREFIX+'/term_profile_'+t+'_tube.pdf')
	plot_terminator_profile(flask_profiles, trace_region, part_region, OUTPUT_PREFIX+'/term_profile_'+t+'_flask.pdf')

########### RIBOZYMES ###########
ribozymes = ['BydvJ', 'PlmJ', 'SarJ', 'RiboJ10', 'RiboJ53', 'RiboJ']
for r in ribozymes:
	trace_region = ['0x58v50', (gff['0x58v50'][r][2]+int(gff['0x58v50'][r][4]['cut_site']))-50, (gff['0x58v50'][r][2]+int(gff['0x58v50'][r][4]['cut_site']))+50]
	part_region = ['0x58v50', gff['0x58v50'][r][2], gff['0x58v50'][r][3]]
	plot_ribozyme_profile(tube_profiles, trace_region, part_region, gff['0x58v50'][r][2]+int(gff['0x58v50'][r][4]['cut_site'])-1, OUTPUT_PREFIX+'/ribo_profile_'+r+'_tube.pdf')
	plot_ribozyme_profile(flask_profiles, trace_region, part_region, gff['0x58v50'][r][2]+int(gff['0x58v50'][r][4]['cut_site'])-1, OUTPUT_PREFIX+'/ribo_profile_'+r+'_flask.pdf')

########### PROMOTERS ###########
promoters = [['pTac', 'pTet1', 'L3S2P55'], ['pBAD1', 'pTet2', 'L3S2P24'], ['pBAD2', 'pBAD2', 'L3S2P11'], ['pBM3R1', 'pAmtR', 'ECK120029600'], ['pSrpR', 'pLitR', 'ECK120033737'], ['pPhlF', 'pPhlF', 'L3S2P21']]
for p in promoters:
	trace_region = ['0x58v50', gff['0x58v50'][p[0]][2]-100, gff['0x58v50'][p[2]][3]]
	part_region = ['0x58v50', gff['0x58v50'][p[0]][2], gff['0x58v50'][p[1]][3]]
	plot_promoter_profile(tube_profiles, trace_region, part_region, OUTPUT_PREFIX+'/pro_profile_'+p[0]+'-'+p[1]+'_tube.pdf')
	plot_promoter_profile(flask_profiles, trace_region, part_region, OUTPUT_PREFIX+'/pro_profile_'+p[0]+'-'+p[1]+'_flask.pdf')
