#!/usr/bin/env python
"""
	Plot transcription profiles for nif transfer experiments
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
import dnaplotlib as dpl

matplotlib.rcParams['lines.dash_joinstyle']  = 'miter'
matplotlib.rcParams['lines.dash_capstyle']   = 'butt'
matplotlib.rcParams['lines.solid_joinstyle'] = 'miter'
matplotlib.rcParams['lines.solid_capstyle']  = 'projecting'
matplotlib.rcParams['pdf.fonttype']          = 42

OUTPUT_PREFIX = './plots/'

plt.rcParams['ytick.major.pad']='5'
plt.rcParams['xtick.major.pad']='5'
fmt_label_size = 8
fmt_edge_width = 2.5
profile_linewidth = 0.6

# Colour maps to use for the genes
col_yellow      = (0.98,0.98,0.227)
col_dark_blue   = (0.129,0.443,0.710)
col_light_blue  = (0.0,0.627,0.878)
col_grey        = (0.839,0.827,0.827)
col_dark_green  = (0.0,0.392,0.196)
col_mid_green   = (0.584,0.749,0.180)
col_light_green = (0.243,0.659,0.227)
col_red         = (0.890,0.184,0.200)
col_black       = (0,0,0)

cmap = {}
cmap['J'] = col_yellow
cmap['H'] = col_dark_blue
cmap['D'] = col_dark_blue
cmap['K'] = col_dark_blue
cmap['T'] = col_grey
cmap['Y'] = col_light_blue
cmap['E'] = col_dark_green
cmap['N'] = col_dark_green
cmap['X'] = col_dark_green
cmap['U'] = col_mid_green
cmap['S'] = col_mid_green
cmap['V'] = col_mid_green
cmap['W'] = col_mid_green
cmap['Z'] = col_mid_green
cmap['M'] = col_mid_green
cmap['F'] = col_yellow
cmap['L'] = col_red
cmap['A'] = col_red
cmap['B'] = col_light_green
cmap['Q'] = col_light_green

cmap['Pwt'] = col_black
cmap['P2'] = col_black
cmap['P2_1'] = col_black
cmap['P2_2'] = col_black
cmap['P3.1'] = col_black
cmap['P4.1'] = col_black
cmap['P4.2'] = col_black

col_sense = (0.5,0.5,0.5)
col_antisense = (0.5,0.5,0.5)

col_sense_low = (0.7,0.7,0.7,0.5)
col_antisense_low = (0.7,0.7,0.7,0.5)
col_sense_high = (0.5,0.5,0.5)
col_antisense_high = (0.5,0.5,0.5)

###############################################################################
# Helper functions
###############################################################################

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
				region = [np.array(new_fwd_profile[new_start_bp:new_end_bp]), 
				          np.array(new_rev_profile[new_start_bp:new_end_bp])]
				break
			else:
				if (start_bp > profile[0] and start_bp < profile[1]) and (end_bp > profile[0] and end_bp < profile[1]):
					fwd_profile = list(profile[2])
					rev_profile = list(profile[3])
					profile_len = len(fwd_profile)
					new_start_bp = start_bp-profile[0]
					new_end_bp = end_bp-profile[0]
					region = [np.array(fwd_profile[new_start_bp:new_end_bp]), 
					          np.array(rev_profile[new_start_bp:new_end_bp])]
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

def load_total_mapped_reads (filename):
	factors = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == 2:
			factors[row[0]] = float(row[2])
	return factors

def load_profile_set (samples, file_prefix):
	profiles = {}
	for s in samples:
		fwd_filename = file_prefix + s + '/' + s + '.fwd.profiles.txt'
		rev_filename = file_prefix + s + '/' + s + '.rev.profiles.txt'
		profiles[s] = load_profiles(fwd_filename, rev_filename)
	return profiles

def gff_to_dnaplotlib (gff, chrom, refactored=False):
	design = []
	for part_name in gff[chrom].keys():
		part_data = gff[chrom][part_name]
		p_name = part_name
		p_type = None
		p_fwd = True
		p_opts = {}
		if part_data[0] == 'gene':
			p_type = 'CDS'
			p_opts['arrowhead_length'] = 100
			p_opts['label'] = part_name
			p_opts['label_y_offset'] = -5.5
			p_opts['label_size'] = 8
			if refactored == True:
				p_opts['hatch'] = '////'
			#p_opts['label_style'] = 'italic'
		if part_data[0] == 'promoter':
			p_type = 'Promoter'
			p_opts['x_extent'] = 160
			p_opts['y_extent'] = 3.5
			p_opts['arrowhead_length'] = 50
			p_opts['linewidth'] = 1.0
		if part_data[0] == 'terminator':
			p_type = 'Terminator'
			p_opts['x_extent'] = 1
			p_opts['y_extent'] = 0.0
			p_opts['linewidth'] = 1.0
		if p_type != None:
			p_start_bp = int(part_data[2])
			p_end_bp = int(part_data[3])
			if part_data[1] != '+':
				p_fwd = False
				p_start_bp = int(part_data[3])
				p_end_bp = int(part_data[2])
			if p_name in cmap.keys():
				p_opts['color'] = cmap[p_name]
			new_part = {'type':p_type,
				        'name':p_name,
				        'start':p_start_bp,
				        'end':p_end_bp,
				        'fwd':p_fwd,
				        'opts':p_opts}
			design.append(new_part)
	return sorted(design, key=lambda k: k['start']) 

def load_plotting_gffs (samples, file_prefix):
	gffs = {}
	for s in samples:
		gff_filename = file_prefix + 'plot_' + s + '.gff'
		gffs[s] = load_gff(gff_filename)
	return gffs

###############################################################################
# Load all the profiles for each organism (specified in BED files)
###############################################################################

#ecoli_samples = ['N1','N2', 'EcoliMG1655_LBWB_1', 'EcoliMG1655_LBWB_2', 
#                 'EcoliMG1655_synnifI4_1', 'EcoliMG1655_synnifI4_2']
ecoli_samples = ['N2', 'EcoliMG1655_LBWB_2', 'EcoliMG1655_synnifI4_1']            
ecoli_file_prefix = '../nif_Ecoli_MG1655/results/'
ecoli_gff_prefix = '../nif_Ecoli_MG1655/data/gff/'
ecoli_profiles = load_profile_set(ecoli_samples, ecoli_file_prefix)
ecoli_gffs = load_plotting_gffs(ecoli_samples, ecoli_gff_prefix)

#rizo_samples = ['N5','N13','N14','N15','N16','N17', 'Rhizobium_1', 'Rhizobium_2', 
#                 'Rhizobium_synnifI4_1', 'Rhizobium_synnifI4_2']
rizo_samples = ['N5','N13','N14','N15','N16','N17', 'Rhizobium_1', 'Rhizobium_2', 
                 'Rhizobium_synnifI4_1', 'Rhizobium_synnifI4_2']
rizo_file_prefix = '../nif_IRBG74_NH/results/'
rizo_gff_prefix = '../nif_IRBG74_NH/data/gff/'
rizo_profiles = load_profile_set(rizo_samples, rizo_file_prefix)
rizo_gffs = load_plotting_gffs(rizo_samples, rizo_gff_prefix)

koxy_samples = ['Koxym5a1_1','Koxym5a1_2', 'Koxym5a1_synnifI4_1', 
                'Koxym5a1_synnifI4_2']
koxy_file_prefix = '../nif_Koxytoca/results/'
koxy_gff_prefix = '../nif_Koxytoca/data/gff/'
koxy_profiles = load_profile_set(koxy_samples, koxy_file_prefix)
koxy_gffs = load_plotting_gffs(koxy_samples, koxy_gff_prefix)

pf5_samples = ['N3','N4','N18','N19','N20','N21']
pf5_file_prefix = '../nif_Pf-5/results/'
pf5_gff_prefix = '../nif_Pf-5/data/gff/'
pf5_profiles = load_profile_set(pf5_samples, pf5_file_prefix)
pf5_gffs = load_plotting_gffs(pf5_samples, pf5_gff_prefix)


#############################################################################
# FIG 1. Koxy5a1-WT
#############################################################################

plot_reversed = True
pad = 100
out_filename = 'fig_01_Koxy5a1-WT.png'
chrom = 'AMPJ01000097.1'
start_bp = (170939-pad)
end_bp = (194229+pad)
profile_data = extract_profile_region(koxy_profiles['Koxym5a1_1'], chrom, start_bp, end_bp)
#norm_factors = load_total_mapped_reads(koxy_file_prefix+'/mapped.reads.matrix.txt')['Koxym5a1_1']
dna_design = gff_to_dnaplotlib(koxy_gffs['Koxym5a1_1'], chrom)
max_y_plus = 22000
max_y_minus = 1800

# Create the figure
fig = plt.figure(figsize=(6.0,2.0))
gs = gridspec.GridSpec(2, 1, height_ratios=[1,0.30])

# Plot the DNA
dr = dpl.DNARenderer(scale=2.0, linewidth=0.8)
ax = plt.subplot(gs[1])
start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
if plot_reversed == True:
	ax.set_xlim([end_bp, start_bp])
else:
	ax.set_xlim([start_bp, end_bp])
ax.set_ylim([-8,5])
ax.plot([start_bp-(pad*10.0),end_bp+(pad*10.0)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Plot the profile
ax = plt.subplot(gs[0])
ax.fill_between(range(start_bp, end_bp), profile_data[0], np.zeros(end_bp-start_bp), facecolor=col_sense, color=col_sense, linewidth=1)
ax.fill_between(range(start_bp, end_bp), profile_data[1]*-1.0, np.zeros(end_bp-start_bp), facecolor=col_antisense, color=col_antisense, linewidth=1)
ax.plot([start_bp, end_bp], [0,0], linewidth=0.8, color=(0,0,0))

# Format the axes and save
ax.set_xticks([])
ax.set_yticks([])
if plot_reversed == True:
	ax.set_xlim([end_bp, start_bp])
	ax.set_ylim([max_y_minus, -max_y_plus])
else:
	ax.set_xlim([start_bp, end_bp])
	ax.set_ylim([-max_y_minus, max_y_plus])
ax.set_yscale('symlog', linthreshy=100, linscaley=0.1)
ax.tick_params(axis='y', which='major', labelsize=fmt_label_size, pad=1, length=2, width=0.5)
for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(0.8)
ax.axis('off')
plt.subplots_adjust(hspace=.00, wspace=.00, left=.01, right=.99, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True, dpi=600)
plt.close('all')

#############################################################################
# FIG 2. Comparison of 7017 WT cluster across hosts
#############################################################################

pad = 100
out_filename = 'fig_02_WT7017_comp.png'
ecoli_chrom = 'pKU7017'
ecoli_norm = 23016036.0/1000000.0#23016036.0
ecoli_start_bp = 2140-pad
ecoli_end_bp = 25614+pad
rizo_chrom = 'pBBR7017'
rizo_norm = 1117055.0/1000000.0# 1117055.0
rizo_start_bp = 7413-pad
rizo_end_bp = 30887+pad
pf5_chrom = 'pPf7017'
pf5_norm = 5430681.0/1000000.0#5430681.0
pf5_start_bp = 7205-pad
pf5_end_bp = 30679+pad
ecoli_profile_data = extract_profile_region(ecoli_profiles['N2'], ecoli_chrom, ecoli_start_bp, ecoli_end_bp)
rizo_profile_data = extract_profile_region(rizo_profiles['N5'], rizo_chrom, rizo_start_bp, rizo_end_bp)
pf5_profile_data = extract_profile_region(pf5_profiles['N3'], pf5_chrom, pf5_start_bp, pf5_end_bp)
dna_design = gff_to_dnaplotlib(ecoli_gffs['N2'], ecoli_chrom)
max_y_plus = 40000 #0.045
max_y_minus = 40000 #0.045 #12

# Create the figure
fig = plt.figure(figsize=(6.0,4.5))
gs = gridspec.GridSpec(4, 1, height_ratios=[0.75,1,0.75,0.32])

# Plot the DNA
dr = dpl.DNARenderer(scale=2.0, linewidth=0.8)
ax = plt.subplot(gs[3])
start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
ax.set_xlim([ecoli_start_bp, ecoli_end_bp])
ax.set_ylim([-8,5])
ax.plot([ecoli_start_bp-(pad*10.0),ecoli_end_bp+(pad*10.0)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Plot the profile
ax = plt.subplot(gs[0])
ax.fill_between(range(ecoli_start_bp, ecoli_end_bp), ecoli_profile_data[0]/ecoli_norm, np.zeros(ecoli_end_bp-ecoli_start_bp), facecolor=col_sense, color=col_sense, linewidth=1)
ax.fill_between(range(ecoli_start_bp, ecoli_end_bp), (ecoli_profile_data[1]*-1.0)/ecoli_norm, np.zeros(ecoli_end_bp-ecoli_start_bp), facecolor=col_antisense, color=col_antisense, linewidth=1)
ax.plot([ecoli_start_bp, ecoli_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([ecoli_start_bp, ecoli_end_bp])
ax.set_ylim([-max_y_minus, max_y_plus])
ax.set_yscale('symlog', linthreshy=100, linscaley=0.1)

ax = plt.subplot(gs[1])
ax.fill_between(range(rizo_start_bp, rizo_end_bp), rizo_profile_data[0]/rizo_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_sense, color=col_sense, linewidth=1)
ax.fill_between(range(rizo_start_bp, rizo_end_bp), (rizo_profile_data[1]*-1.0)/rizo_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_antisense, color=col_antisense, linewidth=1)
ax.plot([rizo_start_bp, rizo_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([rizo_start_bp, rizo_end_bp])
ax.set_yscale('symlog', linthreshy=100, linscaley=0.1)
ax.set_ylim([-max_y_minus, max_y_plus])


ax = plt.subplot(gs[2])
ax.fill_between(range(pf5_start_bp, pf5_end_bp), pf5_profile_data[0]/pf5_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_sense, color=col_sense, linewidth=1)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), (pf5_profile_data[1]*-1.0)/pf5_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_antisense, color=col_antisense, linewidth=1)
ax.plot([pf5_start_bp, pf5_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([pf5_start_bp, pf5_end_bp])
ax.set_ylim([-max_y_minus, max_y_plus])
ax.set_yscale('symlog', linthreshy=100, linscaley=0.1)

# Format the axes and save
plt.subplots_adjust(hspace=.00, wspace=.00, left=.01, right=.99, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True, dpi=600)
plt.close('all')

#############################################################################
# FIG 3. Comparison of v1.1 (MB) cluster across hosts
#############################################################################

pad = 100
out_filename = 'fig_03_MB1.1_comp.png'
rizo_chrom = 'pBBG_MB'
rizo_lowT7_norm = 812740.0
rizo_highT7_norm = 29890466.0
rizo_start_bp = 6956-pad
rizo_end_bp = 26502+pad
pf5_chrom = 'pPf_MB'
pf5_lowT7_norm = 5480126.0
pf5_highT7_norm = 1966434.0
pf5_start_bp = 5626-pad
pf5_end_bp = 25172+pad
rizo_profile_data_lowT7 = extract_profile_region(rizo_profiles['N14'], rizo_chrom, rizo_start_bp, rizo_end_bp)
rizo_profile_data_highT7 = extract_profile_region(rizo_profiles['N16'], rizo_chrom, rizo_start_bp, rizo_end_bp)
pf5_profile_data_lowT7 = extract_profile_region(pf5_profiles['N18'], pf5_chrom, pf5_start_bp, pf5_end_bp)
pf5_profile_data_highT7 = extract_profile_region(pf5_profiles['N20'], pf5_chrom, pf5_start_bp, pf5_end_bp)
dna_design = gff_to_dnaplotlib(rizo_gffs['N14'], rizo_chrom, refactored=False)
max_y_plus = 0.040
max_y_minus = 0.007

# Create the figure
fig = plt.figure(figsize=(6.0,3.0))
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 0.38])

# Plot the DNA
dr = dpl.DNARenderer(scale=2.0, linewidth=0.8)
ax = plt.subplot(gs[2])
start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
ax.set_xlim([rizo_start_bp, rizo_end_bp])
ax.set_ylim([-8,5])
ax.plot([rizo_start_bp-(pad*10.0),rizo_end_bp+(pad*10.0)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Plot the profile
ax = plt.subplot(gs[0])
ax.fill_between(range(rizo_start_bp, rizo_end_bp), rizo_profile_data_lowT7[0]/rizo_lowT7_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_sense_low, color=col_sense_low, linewidth=0, zorder=-1)
ax.fill_between(range(rizo_start_bp, rizo_end_bp), (rizo_profile_data_lowT7[1]*-1.0)/rizo_lowT7_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_antisense_low, color=col_antisense_low, linewidth=0, zorder=-1)
ax.fill_between(range(rizo_start_bp, rizo_end_bp), rizo_profile_data_highT7[0]/rizo_highT7_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_sense_high, color=col_sense_high, linewidth=1, zorder=-2)
ax.fill_between(range(rizo_start_bp, rizo_end_bp), (rizo_profile_data_highT7[1]*-1.0)/rizo_highT7_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_antisense_high, color=col_antisense_high, linewidth=1, zorder=-2)
ax.plot([rizo_start_bp, rizo_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([rizo_start_bp, rizo_end_bp])
#ax.set_yscale('symlog', linthreshx=10)
ax.set_ylim([-max_y_minus, max_y_plus])
ax.tick_params(axis='y', which='major', labelsize=fmt_label_size, pad=1, length=2, width=0.5)

ax = plt.subplot(gs[1])
ax.fill_between(range(pf5_start_bp, pf5_end_bp), pf5_profile_data_lowT7[0]/pf5_lowT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_sense_low, color=col_sense_low, linewidth=0, zorder=-1)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), (pf5_profile_data_lowT7[1]*-1.0)/pf5_lowT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_antisense_low, color=col_antisense_low, linewidth=0, zorder=-1)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), pf5_profile_data_highT7[0]/pf5_highT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_sense_high, color=col_sense_high, linewidth=1, zorder=-2)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), (pf5_profile_data_highT7[1]*-1.0)/pf5_highT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_antisense_high, color=col_antisense_high, linewidth=1, zorder=-2)
ax.plot([pf5_start_bp, pf5_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([pf5_start_bp, pf5_end_bp])
#ax.set_yscale('symlog', linthreshx=10)
ax.set_ylim([-max_y_minus/25.0, max_y_plus/25.0])
ax.tick_params(axis='y', which='major', labelsize=fmt_label_size, pad=1, length=2, width=0.5)

plt.subplots_adjust(hspace=.00, wspace=.00, left=.01, right=.99, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True, dpi=600)
plt.close('all')

#############################################################################
# FIG 4. Comparison of v2.1 (i4) cluster across hosts
#############################################################################

pad = 100
out_filename = 'fig_04_nif60-2.1_comp.png'


ecoli_chrom = 'synI4'
ecoli_highT7_norm = 5471713.0
ecoli_start_bp = 7262-pad
ecoli_end_bp = 26910+pad

rizo_chrom = 'pBBRnif60percent_with_par'
rizo_highT7_norm = 5471713.0
rizo_start_bp = 7262-pad
rizo_end_bp = 26910+pad
pf5_chrom = 'pPf_nif60'
pf5_lowT7_norm = 4746901.0
pf5_highT7_norm = 2032892.0
pf5_start_bp = 5633-pad
pf5_end_bp = 25281+pad

ecoli_profile_data_highT7 = extract_profile_region(rizo_profiles['N17'], rizo_chrom, rizo_start_bp, rizo_end_bp)

rizo_profile_data_highT7 = extract_profile_region(rizo_profiles['N17'], rizo_chrom, rizo_start_bp, rizo_end_bp)
pf5_profile_data_lowT7 = extract_profile_region(pf5_profiles['N19'], pf5_chrom, pf5_start_bp, pf5_end_bp)
pf5_profile_data_highT7 = extract_profile_region(pf5_profiles['N21'], pf5_chrom, pf5_start_bp, pf5_end_bp)
dna_design = gff_to_dnaplotlib(rizo_gffs['N17'], rizo_chrom, refactored=False)
max_y_plus = 0.005
max_y_minus = 0.0015

# Create the figure
fig = plt.figure(figsize=(6.0,3.0))
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 0.38])

# Plot the DNA
dr = dpl.DNARenderer(scale=2.0, linewidth=0.8)
ax = plt.subplot(gs[2])
start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
ax.set_xlim([rizo_start_bp, rizo_end_bp])
ax.set_ylim([-8,5])
ax.plot([rizo_start_bp-(pad*10.0),rizo_end_bp+(pad*10.0)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Plot the profile
ax = plt.subplot(gs[0])
ax.fill_between(range(rizo_start_bp, rizo_end_bp), rizo_profile_data_highT7[0]/rizo_highT7_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_sense_high, color=col_sense_high, linewidth=1, zorder=-2)
ax.fill_between(range(rizo_start_bp, rizo_end_bp), (rizo_profile_data_highT7[1]*-1.0)/rizo_highT7_norm, np.zeros(rizo_end_bp-rizo_start_bp), facecolor=col_antisense_high, color=col_antisense_high, linewidth=1, zorder=-2)
ax.plot([rizo_start_bp, rizo_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([rizo_start_bp, rizo_end_bp])
#ax.set_yscale('symlog', linthreshx=10)
ax.set_ylim([-max_y_minus, max_y_plus])
ax.tick_params(axis='y', which='major', labelsize=fmt_label_size, pad=1, length=2, width=0.5)

ax = plt.subplot(gs[1])
ax.fill_between(range(pf5_start_bp, pf5_end_bp), pf5_profile_data_lowT7[0]/pf5_lowT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_sense_low, color=col_sense_low, linewidth=0, zorder=-1)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), (pf5_profile_data_lowT7[1]*-1.0)/pf5_lowT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_antisense_low, color=col_antisense_low, linewidth=0, zorder=-1)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), pf5_profile_data_highT7[0]/pf5_highT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_sense_high, color=col_sense_high, linewidth=1, zorder=-2)
ax.fill_between(range(pf5_start_bp, pf5_end_bp), (pf5_profile_data_highT7[1]*-1.0)/pf5_highT7_norm, np.zeros(pf5_end_bp-pf5_start_bp), facecolor=col_antisense_high, color=col_antisense_high, linewidth=1, zorder=-2)
ax.plot([pf5_start_bp, pf5_end_bp], [0,0], linewidth=1.0, color=(0,0,0))
ax.axis('off')
ax.set_xticks([])
ax.set_xlim([pf5_start_bp, pf5_end_bp])
#ax.set_yscale('symlog', linthreshx=10)
ax.set_ylim([-max_y_minus/1.0, max_y_plus/1.0])
ax.tick_params(axis='y', which='major', labelsize=fmt_label_size, pad=1, length=2, width=0.5)

plt.subplots_adjust(hspace=.00, wspace=.00, left=.01, right=.99, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True, dpi=600)
plt.close('all')



