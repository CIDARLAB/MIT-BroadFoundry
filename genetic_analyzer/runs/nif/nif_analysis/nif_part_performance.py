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

plt.rcParams['ytick.major.pad']='0.5'
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

cmap['P7'] = col_black
cmap['P7.2'] = col_black
cmap['P7.2_1'] = col_black
cmap['P7.2_2'] = col_black
cmap['P7.3.1'] = col_black
cmap['P7.4.1'] = col_black
cmap['P7.4.2'] = col_black

col_sense = (0.5,0.5,0.5)
col_antisense = (0.5,0.5,0.5)

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

def gff_to_dnaplotlib (gff, chrom):
	design = []
	for part_name in gff[chrom].keys():
		part_data = gff[chrom][part_name]
		p_name = part_name
		p_type = None
		p_fwd = True
		p_opts = {}
		if part_data[0] == 'gene':
			p_type = 'CDS'
			#p_opts['arrowhead_length'] = 100
			p_opts['label'] = part_name
			p_opts['label_y_offset'] = -1.0
			p_opts['label_x_offset'] = -2.0
			p_opts['label_size'] = 7
			p_opts['y_extent'] = 7
			p_opts['y_extent'] = 7
			p_opts['arrowhead_height'] = 4
			p_opts['arrowhead_length'] = 10
			p_opts['label_style'] = 'italic'
		if part_data[0] == 'promoter':
			p_type = 'Promoter'
			p_opts['label'] = part_name.split('_')[0]
			p_opts['label_y_offset'] = 20.0
			p_opts['label_size'] = 8
			p_opts['linewidth'] = 1.0
		if part_data[0] == 'terminator':
			p_type = 'Terminator'
			p_opts['label'] = part_name.split('_')[0]
			p_opts['label_y_offset'] = -13.0
			p_opts['label_x_offset'] = 4.0
			p_opts['label_size'] = 8
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
# FIG 5. Nif variant designs
#############################################################################

pad = 100
out_filename = 'fig_05_design_variants.pdf'

dna_design_v1_0 = gff_to_dnaplotlib(rizo_gffs['N13'], 'pBBR_OR')
start_bp_v1_0 = (7251-pad)
end_bp_v1_0 = (25587+pad)
dna_design_v1_1_MB = gff_to_dnaplotlib(rizo_gffs['N14'], 'pBBG_MB')
start_bp_v1_1_MB = (6956-pad)
end_bp_v1_1_MB = (26502+pad)
dna_design_v2_1_i4 = gff_to_dnaplotlib(rizo_gffs['N17'], 'pBBRnif60percent_with_par')
start_bp_v2_1_i4 = (7262-pad)
end_bp_v2_1_i4 = (26910+pad)

# Create the figure
fig = plt.figure(figsize=(4.4,1.5))
gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,1])
min_x = -32
max_x = 32

# Plot v1.0
dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
ax = plt.subplot(gs[0])
start, end = dr.renderDNA(ax, dna_design_v1_0, dr.SBOL_part_renderers())
ax.set_xlim([start, end])
ax.set_ylim([min_x,max_x])
ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Plot v1.1 MB
dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
ax = plt.subplot(gs[1])
start, end = dr.renderDNA(ax, dna_design_v1_1_MB, dr.SBOL_part_renderers())
ax.set_xlim([start, end])
ax.set_ylim([min_x,max_x])
ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Plot v2.1 i4
dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
ax = plt.subplot(gs[2])
start, end = dr.renderDNA(ax, dna_design_v2_1_i4, dr.SBOL_part_renderers())
ax.set_xlim([start, end])
ax.set_ylim([min_x,max_x])
ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Save the figure
plt.subplots_adjust(hspace=.00, wspace=.00, left=.05, right=.95, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True)
plt.close('all')

#############################################################################
# FIG 6. Terminator performance
#############################################################################

def add_terminator_profile (ax, profiles, profile_cols):
	fmt_char_data_linewidth = 1.0
	fmt_axis_outline_width = 1.2
	fmt_char_line_width = 3.0
	fmt_label_size = 9
	annot_line_col = (0,0,0)
	for s_idx in range(len(profiles)):
		ax.plot(range(len(profiles[s_idx])), profiles[s_idx], color=profile_cols[s_idx], alpha=1, linewidth=1.5)
		#ax.plot(range(len(profiles[s_idx])), profiles[s_idx], color=(0,0,0), alpha=1, linewidth=1.5)
	plt.axvspan(20, len(profiles[0])-20, facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([0,len(profiles[0])])
	ax.set_xscale('linear')
	ax.set_ylim([0,100000])
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

def get_part_profile (profiles, gff, sample, chrom, part_name, us_pad=20, ds_pad=20):
	start_bp = gff[sample][chrom][part_name][2]-us_pad
	end_bp = gff[sample][chrom][part_name][3]+ds_pad
	return extract_profile_region(profiles[sample], chrom, start_bp, end_bp)

fig = plt.figure(figsize=(6.5,1.2))
gs = gridspec.GridSpec(1, 5)
pcols = [col_light_blue, col_light_blue, col_mid_green, col_mid_green]
ax = plt.subplot(gs[0])
Pf5_terms = []
Pf5_terms.append(get_part_profile(pf5_profiles, pf5_gffs, 'N19', 'pPf_nif60', 'T7')[0])
Pf5_terms.append(get_part_profile(pf5_profiles, pf5_gffs, 'N21', 'pPf_nif60', 'T7')[0])
add_terminator_profile(ax, Pf5_terms, pcols)
for t_idx, t in enumerate(['T4', 'T6_1', 'T6_2', 'Twt']):
	ax = plt.subplot(gs[t_idx+1])	
	Pf5_terms = []
	Pf5_terms.append(get_part_profile(pf5_profiles, pf5_gffs, 'N19', 'pPf_nif60', t)[0])
	Pf5_terms.append(get_part_profile(pf5_profiles, pf5_gffs, 'N21', 'pPf_nif60', t)[0])
	Pf5_terms.append(get_part_profile(pf5_profiles, pf5_gffs, 'N18', 'pPf_MB', t)[0])
	Pf5_terms.append(get_part_profile(pf5_profiles, pf5_gffs, 'N20', 'pPf_MB', t)[0])
	add_terminator_profile(ax, Pf5_terms, pcols)
plt.subplots_adjust(hspace=.00, wspace=.32, left=.05, right=.99, top=.90, bottom=.10)
plt.savefig(OUTPUT_PREFIX+'fig_06_term_traces_Pf-5.pdf', transparent=True)
plt.close('all')

fig = plt.figure(figsize=(6.5,1.2))
gs = gridspec.GridSpec(1, 5)
pcols = [col_light_blue, col_light_blue, col_mid_green, col_mid_green, col_red]
ax = plt.subplot(gs[0])
rizo_terms = []
rizo_terms.append(get_part_profile(rizo_profiles, rizo_gffs, 'N17', 'pBBRnif60percent_with_par', 'T7')[0])
add_terminator_profile(ax, rizo_terms, pcols)
for t_idx, t in enumerate(['T4', 'T6_1', 'T6_2', 'Twt']):
	ax = plt.subplot(gs[t_idx+1])	
	rizo_terms = []
	rizo_terms.append(get_part_profile(rizo_profiles, rizo_gffs, 'N17', 'pBBRnif60percent_with_par', t)[0])
	rizo_terms.append(get_part_profile(rizo_profiles, rizo_gffs, 'N13', 'pBBR_OR', t)[0])
	rizo_terms.append(get_part_profile(rizo_profiles, rizo_gffs, 'N15', 'pBBR_OR', t)[0])
	rizo_terms.append(get_part_profile(rizo_profiles, rizo_gffs, 'N14', 'pBBG_MB', t)[0])
	rizo_terms.append(get_part_profile(rizo_profiles, rizo_gffs, 'N16', 'pBBG_MB', t)[0])
	add_terminator_profile(ax, rizo_terms, pcols)
plt.subplots_adjust(hspace=.00, wspace=.32, left=.05, right=.99, top=.90, bottom=.10)
plt.savefig(OUTPUT_PREFIX+'fig_06_term_traces_rizo.pdf', transparent=True)
plt.close('all')

#############################################################################
# FIG 7. Promoter performance
#############################################################################

def add_promoter_profile (ax, profiles, profile_cols):
	fmt_char_data_linewidth = 1.0
	fmt_axis_outline_width = 1.2
	fmt_char_line_width = 3.0
	fmt_label_size = 9
	annot_line_col = (0,0,0)
	for s_idx in range(len(profiles)):
		ax.plot(range(len(profiles[s_idx])), profiles[s_idx], color=profile_cols[s_idx], alpha=1, linewidth=1.5)
		#ax.plot(range(len(profiles[s_idx])), profiles[s_idx], color=(0,0,0), alpha=1, linewidth=1.5)
	plt.axvspan(20, len(profiles[0])-20, facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([0,len(profiles[0])])
	ax.set_xscale('linear')
	ax.set_ylim([0,100000])
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

def load_norm_factors (filename):
	data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == 3:
			chrom = row[0]
			norm_fac = float(row[1])
			lib_size = float(row[2])
			data[chrom] = [lib_size*norm_fac, norm_fac, lib_size]
	return data

Pf5_norm_facs = load_norm_factors('../nif_IRBG74_NH/results/')



