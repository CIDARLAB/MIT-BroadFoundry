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
from matplotlib.ticker import MultipleLocator
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

tube_samples = ['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8']
flask_samples = ['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8']

plt.rcParams['ytick.major.pad']='1' # 5 for all but repressor graphs
plt.rcParams['xtick.major.pad']='5'
fmt_label_size = 13.5
fmt_edge_width = 2.5
trace_height = 2.5

col_sense_rd = (0.5,0.5,0.5)
col_antisense_rd = (0.9,0.5,0.5)

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
	fig = plt.figure(figsize=(3.0,8.5))
	gs = gridspec.GridSpec(8, 1, height_ratios=[1,1,1,1,1,1,1,1])
	for s in profiles.keys():
		s_idx = int(s.split('_')[1])-1
		ax = plt.subplot(gs[s_idx])
		profile_data = extract_profile_region(profiles[s], trace_region[0], trace_region[1], trace_region[2])
		plt.axvspan(part_region[1], part_region[2], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
		ax.fill_between(range(trace_region[1], trace_region[2]), np.array(profile_data[0]), np.zeros(len(profile_data[0])), facecolor=col_sense_rd, linewidth=0)
		ax.fill_between(range(trace_region[1], trace_region[2]), np.array(profile_data[1])*-1.0, np.zeros(len(profile_data[1])), facecolor=col_antisense_rd, linewidth=0)
		ax.plot([trace_region[1], trace_region[2]], [0,0], color=(0.0,0.0,0.0), alpha=1, linewidth=1.5)
		ax.tick_params(axis='x', labelsize=fmt_label_size)
		ax.tick_params(axis='y', labelsize=fmt_label_size)
		ax.set_xlim([trace_region[1],trace_region[2]])
		ax.set_xscale('linear')
		ax.set_ylim([-100000,100000])
		ax.set_yscale('symlog', linthreshy=10)
		ax.yaxis.set_ticklabels([])
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.tick_params(axis='x', labelsize=fmt_label_size)
		ax.tick_params(axis='y', labelsize=fmt_label_size)
		ax.set_xticks([])
		ax.set_xticklabels([], visible=False)
	plt.subplots_adjust(hspace=0.0, left=0.01, right=0.99, top=0.99, bottom=0.01)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')

def plot_antisense_stats (profiles, trace_regions, part_regions, filename_out, point_shape='o'):
	us_region_bp = 200
	ds_region_bp = 50
	fracs_antisense = []
	term_strengths = []
	
	for part_idx in range(len(part_regions)):
		trace_region = trace_regions[part_idx]
		part_region = part_regions[part_idx]
		for s in profiles.keys():
			if s not in ['L3S3P21-2', 'L3S2P22']:
				profile_data = extract_profile_region(profiles[s], trace_region[0], trace_region[1], trace_region[2])
				if float(np.mean(profile_data[0][0:us_region_bp])) == 0.0:
					fracs_antisense.append(0.0)
				else:
					fracs_antisense.append(float(np.mean(profile_data[1][us_region_bp:-ds_region_bp]))/float(np.mean(profile_data[0][0:us_region_bp])))
				# Calculate termination strength statistics
				us_reads = np.mean(profile_data[0][0:us_region_bp])
				ds_reads = np.mean(profile_data[0][-ds_region_bp:])
				t_s = 0.0
				t_e = 0.0
				t_type = 1
				if ds_reads == 0 and us_reads > 0:
					t_type = 2
					t_e = 1.0/us_reads
				else:
					t_e = ds_reads/us_reads				
				if us_reads < 250:
					t_type = 0
				if t_e != 0.0:
					t_s = 1.0/t_e
				term_strengths.append([t_s, t_type])

	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.8,2.4))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])

	for idx in range(len(fracs_antisense)):
		x_val = term_strengths[idx][0]
		y_val = fracs_antisense[idx]
		t_type = term_strengths[idx][1]
		point_col = (0,0,0)
		p_zorder = 10
		if t_type == 0:
			#point_col = (0.7,0.7,0.7)
			p_zorder = 9
		#if t_type == 2:
			#point_col = (0.95,0.30,0.25)
		ax.scatter([x_val], [y_val], marker=point_shape, zorder=p_zorder, s=40, color=(1,1,1,0), edgecolor=point_col, linewidth=2.0)
		print x_val, y_val

	ax.set_xlim([1, 10000])
	ax.set_ylim([-0.03, 50])
	ax.set_xscale('symlog', linthreshy=0.1)
	ax.set_yscale('symlog', linthreshy=0.1)

	plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.15)
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
trace_regions = []
part_regions = []
for t in terminators:
	trace_region = ['0x58v50', gff['0x58v50'][t][2]-100, gff['0x58v50'][t][3]+100]
	trace_region_stats = ['0x58v50', gff['0x58v50'][t][2]-200, gff['0x58v50'][t][3]+50]
	part_region = ['0x58v50', gff['0x58v50'][t][2], gff['0x58v50'][t][3]]
	plot_terminator_profile(tube_profiles, trace_region, part_region, OUTPUT_PREFIX+'/antisense_term_profile_'+t+'_tube.pdf')
	plot_terminator_profile(flask_profiles, trace_region, part_region, OUTPUT_PREFIX+'/antisense_term_profile_'+t+'_flask.pdf')
	trace_regions.append(trace_region_stats)
	part_regions.append(part_region)

plot_antisense_stats(tube_profiles, trace_regions, part_regions, OUTPUT_PREFIX+'/antisense_term_stats_tube.pdf', point_shape='o')
plot_antisense_stats(flask_profiles, trace_regions, part_regions, OUTPUT_PREFIX+'/antisense_term_stats_flask.pdf', point_shape='x')
