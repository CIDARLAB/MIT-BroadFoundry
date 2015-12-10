#!/usr/bin/env python
"""
	Plot transcription profile plots for paper.
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import dnaplotlib as dpl

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

# Colour maps to use for the genes
cmap = {}
cmap['AmtR'] = (1.0,0.75,0.17) # 255, 193, 43
cmap['LitR'] = (0.38,0.82,0.32) # 98, 209, 83
cmap['BM3R1'] = (0.95,0.30,0.25) # 242, 78, 65
cmap['SrpR'] = (0.38,0.65,0.87) # 97, 165, 223
cmap['PhlF'] = (0.55,0.35,0.64) # 141, 89, 163
cmap['YFP'] = (0.98,0.97,0.35) # 250, 248, 89

cmap['pAmtR'] = (1.0,0.75,0.17) # 255, 193, 43
cmap['pLitR'] = (0.38,0.82,0.32) # 98, 209, 83
cmap['pBM3R1'] = (0.95,0.30,0.25) # 242, 78, 65
cmap['pSrpR'] = (0.38,0.65,0.87) # 97, 165, 223
cmap['pPhlF'] = (0.55,0.35,0.64) # 141, 89, 163
cmap['pYFP'] = (0.98,0.97,0.35) # 250, 248, 89

cmap['pTac'] = (0,0,0)
cmap['pTet1'] = (0,0,0)
cmap['pBAD1'] = (0,0,0)
cmap['pTet2'] = (0,0,0)
cmap['pBAD2'] = (0,0,0)

# Terminator colours (all black)
cmap['L3S2P55'] = (0.0,0.0,0.0) # AmtR
cmap['L3S2P24'] = (0.0,0.0,0.0) # LitR
cmap['L3S2P11'] = (0.0,0.0,0.0) # BM3R1
cmap['ECK120029600'] = (0.0,0.0,0.0) # SrpR
cmap['ECK120033737'] = (0.0,0.0,0.0) # PhlF
cmap['L3S2P21'] = (0.0,0.0,0.0) # YFP

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

def load_norm_factors (filename):
	factors = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == 3:
			factors[row[0]] = float(row[1])*float(row[2])
	return factors

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
		ax.set_yscale('symlog', linthreshy=50, linscaley=0.1)
		ax.yaxis.set_ticklabels([])
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.tick_params(axis='x', labelsize=fmt_label_size)
		ax.tick_params(axis='y', labelsize=fmt_label_size)
		if s_idx == 7:
			print 'here'
		else:
			ax.set_xticks([])
			ax.set_xticklabels([], visible=False)
	plt.subplots_adjust(hspace=0.0, left=0.1, right=0.99, top=0.99, bottom=0.05)
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
			#if s not in ['L3S3P21-2', 'L3S2P22']:
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

	fmt_char_data_linewidth = 1.0
	fmt_axis_outline_width = 1.0
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(1.75,1.6))
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
		ax.scatter([x_val], [y_val], marker=point_shape, zorder=p_zorder, s=15, color=(1,1,1,0), edgecolor=point_col, linewidth=1.2)
		print x_val, y_val

	ax.set_xlim([1, 10000])
	ax.set_ylim([-0.03, 50])
	ax.set_xscale('symlog', linthreshy=0.1)
	ax.set_yscale('symlog', linthreshy=0.1)
	ax.tick_params(axis='y', which='major', labelsize=8, pad=1, length=2, width=0.5)
	ax.tick_params(axis='x', which='major', labelsize=8, pad=3, length=2, width=0.5)
	for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)

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

########### Read-through ###########

def plot_read_through (profiles, pre_gene_ranges, filename, point_shape='o'):
	fmt_char_data_linewidth = 1.0
	fmt_axis_outline_width = 1.0
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(1.75,1.3))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])
	
	# Calculate the read-through
	
	x_vals = [[],[],[],[],[],[]]
	y_vals = [[],[],[],[],[],[]]

	for p_key in profiles.keys():
		p = np.array(profiles[p_key][0])
		for g_idx in range(len(pre_gene_ranges)):
			g_range = pre_gene_ranges[g_idx]
			# Calculate average read-through level
			x_vals[g_idx].append(g_idx)
			y_vals[g_idx].append(np.mean(p[g_range[0]:g_range[1]]))
	for g_idx in range(len(x_vals)):
		ax.scatter(x_vals[g_idx], y_vals[g_idx], marker=point_shape, zorder=-1, s=15, color=(1,1,1,0), edgecolor=(0,0,0), linewidth=1.2)
		mean_rt = np.mean(y_vals[g_idx])
		ax.plot([g_idx-0.32, g_idx+0.32], [mean_rt, mean_rt], color=(0.5,0.5,0.5), linewidth=2.25, zorder=-100)

	plt.axvline(x=0.5, linewidth=0.6, color=(0,0,0))
	plt.axvline(x=1.5, linewidth=0.6, color=(0,0,0))
	plt.axvline(x=2.5, linewidth=0.6, color=(0,0,0))
	plt.axvline(x=3.5, linewidth=0.6, color=(0,0,0))
	plt.axvline(x=4.5, linewidth=0.6, color=(0,0,0))

	ax.set_xlim([-0.5, 5.5])
	ax.set_ylim([-5, 5000])
	ax.set_yscale('symlog', linthreshy=10)
	ax.tick_params(axis='y', which='major', labelsize=8, pad=1, length=2, width=0.5)
	ax.tick_params(axis='x', which='major', labelsize=8, pad=3, length=2, width=0.5)
	plt.setp(ax.get_xticklabels(), visible=False)
	for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)

	plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.15)
	plt.savefig(filename, transparent=True)
	plt.close('all')

########### Predicted promoter sites ###########

def gff_to_dnaplotlib (gff, chrom):
	design = []
	for part_name in gff[chrom].keys():
		part_data = gff[chrom][part_name]
		p_name = part_name
		p_type = None
		p_opts = {}
		if part_data[0] == 'gene':
			p_type = 'CDS'
			p_opts['label'] = part_name
			p_opts['label_y_offset'] = -6.5
			p_opts['label_size'] = 8
			p_opts['label_style'] = 'italic'
		if part_data[0] == 'promoter':
			p_type = 'Promoter'
		if part_data[0] == 'terminator':
			p_type = 'Terminator'
		if p_type != None:
			p_start_bp = int(part_data[2])
			p_end_bp = int(part_data[3])
			if p_name in cmap.keys():
				p_opts['color'] = cmap[p_name]
			new_part = {'type':p_type,
				        'name':p_name,
				        'start':p_start_bp,
				        'end':p_end_bp,
				        'fwd':True,
				        'opts':p_opts}
			design.append(new_part)
	return sorted(design, key=lambda k: k['start']) 

def load_sigma70_sites (filename, plasmid_len=11115, offset=19):
	sites = []
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process each line
	for row in data_reader:
		if len(row) == 20 and row[0] != '# Sequence':
			# Extract the details
			site_dir = '+'
			if '_r' in row[0]:
				site_dir = '-'
			tss = 0
			if site_dir == '+':
				tss = int(row[8])+7+int(row[15])-offset
			else:
				tss = plasmid_len-(int(row[8])-7-int(row[15]))-offset
			tot_score = float(row[17])
			z_score = float(row[18])
			sites.append([site_dir, tss, tot_score, z_score])
	return sites	

def plot_profile_with_sites (profile, sites, terminator_ranges, dna_design, filename, y_range=[-1000000, 1000000]):
	fmt_char_data_linewidth = 1.0
	fmt_axis_outline_width = 1.0
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(3.5,1.0))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1,0.6])
	
	# Plot the DNA
	dr = dpl.DNARenderer(scale=2, linewidth=0.75)
	ax = plt.subplot(gs[1])
	start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
	ax.set_xlim([0, 6124])
	ax.set_ylim([-10,8])
	ax.plot([0,6124], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')
	# Plot the trace
	ax = plt.subplot(gs[0])
	ax.fill_between(range(len(profile[0])), profile[0], np.zeros(len(profile[0])), facecolor=col_sense_rd, linewidth=0)
	ax.fill_between(range(len(profile[1])), np.array(profile[1])*-1.0, np.zeros(len(profile[1])), facecolor=col_antisense_rd, linewidth=0)
	ax.plot([0,len(profile[1])],[0,0], color=(0,0,0), linewidth=1.0)
	# Annotate the sigma-70 sites
	for s in sites:
		if s[0] == '-' and s[3] >= 0.5: #s[2] > 0.0 and 
			ax.annotate(r'$\blacktriangle$',(s[1], 0), ha='center', va='top', fontsize=6)
	# Terminator ranges
	for t_range in terminator_ranges:
		ax.fill_between([t_range[0], t_range[1]], [-9999999999, -9999999999], [9999999999, 9999999999], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
	ax.set_xlim([0, len(profile[0])])
	ax.set_xticks([])
	ax.set_ylim(y_range)
	plt.setp(ax.get_yticklabels(), visible=False)
	ax.set_yscale('symlog', linthreshy=50, linscaley=0.5)
	ax.tick_params(axis='y', which='major', labelsize=8, pad=1, length=2, width=0.5)
	for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(hspace=0.01, left=0.01, right=0.99, top=0.96, bottom=0.04)
	plt.savefig(filename, transparent=True)
	plt.close('all')

def plot_profile_set_with_sites (profiles, sites, terminator_ranges, dna_design, filename, y_range=[-1000000, 1000000]):
	fmt_char_data_linewidth = 1.0
	fmt_axis_outline_width = 1.0
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(4.0,5.0))
	gs = gridspec.GridSpec(9, 1, height_ratios=[1,1,1,1,1,1,1,1,0.6])
	
	# Plot the DNA
	dr = dpl.DNARenderer(scale=2, linewidth=0.75)
	ax = plt.subplot(gs[8])
	start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
	ax.set_xlim([0, 6124])
	ax.set_ylim([-10,8])
	ax.plot([0,6124], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')
	
	# Plot the traces
	ax_idx = 0
	for p in sorted(profiles.keys()):
		ax = plt.subplot(gs[ax_idx])
		ax_idx += 1
		ax.fill_between(range(len(profiles[p][0])), profiles[p][0], np.zeros(len(profiles[p][0])), facecolor=col_sense_rd, linewidth=0)
		ax.fill_between(range(len(profiles[p][1])), np.array(profiles[p][1])*-1.0, np.zeros(len(profiles[p][1])), facecolor=col_antisense_rd, linewidth=0)
		ax.plot([0,len(profiles[p][1])],[0,0], color=(0,0,0), linewidth=1.0)
		# Annotate the sigma-70 sites
		for s in sites:
			if s[0] == '-' and s[3] >= 0.5: #s[2] > 0.0 and 
				ax.annotate(r'$\blacktriangle$',(s[1], 0), ha='center', va='top', fontsize=6)
		# Terminator ranges
		for t_range in terminator_ranges:
			ax.fill_between([t_range[0], t_range[1]], [-9999999999, -9999999999], [9999999999, 9999999999], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
		ax.set_xlim([0, len(profiles[p][0])])
		ax.set_xticks([])
		ax.set_ylim(y_range)
		plt.setp(ax.get_yticklabels(), visible=False)
		ax.set_yscale('symlog', linthreshy=50, linscaley=0.5)
		ax.tick_params(axis='y', which='major', labelsize=8, pad=1, length=2, width=0.5)
		for axis in ['top','bottom','left','right']:
				ax.spines[axis].set_linewidth(fmt_axis_outline_width)
				ax.spines[axis].set_linewidth(fmt_axis_outline_width)
				ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	
	plt.subplots_adjust(hspace=0.0, left=0.01, right=0.99, top=0.99, bottom=0.01)
	plt.savefig(filename, transparent=True)
	plt.close('all')

# Calculate a conversion factor from FPKM to reads use tube_4 as high number of reads and YFP FPKM
norm_factors = load_norm_factors(RESULTS_PREFIX+'/norm.factors.matrix.txt')
YFP_reads = np.array(extract_profile_region(tube_profiles['tube_4'], '0x58v50', gff['0x58v50']['YFP'][2]+100, gff['0x58v50']['YFP'][3]-100)[0])/norm_factors['tube_4']
YFP_FPKM = 38977.39525
reads_FPKM_fac = YFP_FPKM/np.median(YFP_reads)
dna_design = gff_to_dnaplotlib(gff, '0x58v50')

terminator_ranges = []
for t in ['ECK120029600', 'ECK120033737', 'L3S2P11', 'L3S2P21', 'L3S2P24', 'L3S2P55']:
	terminator_ranges.append([gff['0x58v50'][t][2], gff['0x58v50'][t][3]])

# Load all the predicted sigma-70 sites for circuit
sig70_sites = load_sigma70_sites('./sigma70_scanner/0x58v50_results.txt', plasmid_len=11115, offset=19)
cir_tube_profile_data = {}
for s in tube_samples:
	extracted_profile = extract_profile_region(tube_profiles[s], '0x58v50', gff['0x58v50']['pTac'][2], gff['0x58v50']['L3S2P21'][3])
	new_profile = []
	new_profile.append((np.array(extracted_profile[0])/norm_factors[s])*reads_FPKM_fac)
	new_profile.append((np.array(extracted_profile[1])/norm_factors[s])*reads_FPKM_fac)
	cir_tube_profile_data[s] = new_profile
	
cir_flask_profile_data = {}
for s in flask_samples:
	extracted_profile = extract_profile_region(flask_profiles[s], '0x58v50', gff['0x58v50']['pTac'][2], gff['0x58v50']['L3S2P21'][3]) 
	new_profile = []
	new_profile.append((np.array(extracted_profile[0])/norm_factors[s])*reads_FPKM_fac)
	new_profile.append((np.array(extracted_profile[1])/norm_factors[s])*reads_FPKM_fac)
	cir_flask_profile_data[s] = new_profile

pre_gene_ranges = []
pre_gene_ranges.append([0, gff['0x58v50']['pTac'][2]])
pre_gene_ranges.append([gff['0x58v50']['L3S2P55'][3], gff['0x58v50']['pBAD1'][2]])
pre_gene_ranges.append([gff['0x58v50']['L3S2P24'][3], gff['0x58v50']['pBAD2'][2]])
pre_gene_ranges.append([gff['0x58v50']['L3S2P11'][3], gff['0x58v50']['pBM3R1'][2]])
pre_gene_ranges.append([gff['0x58v50']['ECK120029600'][3], gff['0x58v50']['pSrpR'][2]])
pre_gene_ranges.append([gff['0x58v50']['ECK120033737'][3], gff['0x58v50']['pPhlF'][2]])
plot_read_through(cir_tube_profile_data, pre_gene_ranges, './plots_general/tube_read_through.pdf')

# Plot for Figure 5
plot_profile_with_sites(cir_tube_profile_data['tube_1'], sig70_sites, terminator_ranges, dna_design, './plots_general/fig_4_tube_1.pdf', y_range=[-100000, 100000])

# Plot SI figure
plot_profile_set_with_sites(cir_tube_profile_data, sig70_sites, terminator_ranges, dna_design, './plots_general/tube_antisense.pdf', y_range=[-1000000, 1000000])
plot_profile_set_with_sites(cir_flask_profile_data, sig70_sites, terminator_ranges, dna_design, './plots_general/flask_antisense.pdf', y_range=[-1000000, 1000000])


########### TERMINATORS ###########
terminators = ['ECK120029600', 'ECK120033737', 'L3S2P11', 'L3S2P21', 'L3S2P24', 'L3S2P55', 'L3S2P22', 'L3S3P21-2']
trace_regions = []
part_regions = []
for t in terminators:
	trace_region = ['0x58v50', gff['0x58v50'][t][2]-1000, gff['0x58v50'][t][3]+1000]
	trace_region_stats = ['0x58v50', gff['0x58v50'][t][2]-200, gff['0x58v50'][t][3]+50]
	part_region = ['0x58v50', gff['0x58v50'][t][2], gff['0x58v50'][t][3]]
	plot_terminator_profile(tube_profiles, trace_region, part_region, OUTPUT_PREFIX+'/antisense_term_profile_'+t+'_tube.pdf')
	plot_terminator_profile(flask_profiles, trace_region, part_region, OUTPUT_PREFIX+'/antisense_term_profile_'+t+'_flask.pdf')
	trace_regions.append(trace_region_stats)
	part_regions.append(part_region)

plot_antisense_stats(tube_profiles, trace_regions, part_regions, OUTPUT_PREFIX+'/antisense_term_stats_tube.pdf', point_shape='o')
plot_antisense_stats(flask_profiles, trace_regions, part_regions, OUTPUT_PREFIX+'/antisense_term_stats_flask.pdf', point_shape='x')

########### CALCULATE ANTISENSE FRACTION ###########
print 'Tube antisense frac:'
for s in tube_samples:
	cir_profile_data = extract_profile_region(tube_profiles[s], '0x58v50', gff['0x58v50']['pTac'][2], gff['0x58v50']['L3S2P21'][3])
	print s, np.sum(cir_profile_data[1])/np.sum(cir_profile_data[0])
print 'Flask antisense frac:'
for s in flask_samples:
	cir_profile_data = extract_profile_region(flask_profiles[s], '0x58v50', gff['0x58v50']['pTac'][2], gff['0x58v50']['L3S2P21'][3])
	print s, np.sum(cir_profile_data[1])/np.sum(cir_profile_data[0])

