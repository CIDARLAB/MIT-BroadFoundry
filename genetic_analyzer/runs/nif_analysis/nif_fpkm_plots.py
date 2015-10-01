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

def load_fpkms (filename):
	fpkms = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	header = next(data_reader)
	header_map = {}
	for idx, h in enumerate(header):
		header_map[idx] = header[idx]
		fpkms[header[idx]] = {}
	# Process each line
	for row in data_reader:
		if len(row) > 1:
			tag = row[0]
			vals = [float(x) for x in row[1:]]
			for idx, v in enumerate(vals):
				fpkms[header_map[idx+1]][tag] = v 
	return fpkms

def find_nif_key (fpkms, nif_name):
	for cur_key in fpkms.keys():
		if nif_name in cur_key:
			return cur_key
	return None

def format_dna_design_fpkm (dna_design, fpkms, cmapper, norm):
	for part in dna_design:
		if part['type'] == 'CDS':
			nif_key = find_nif_key (fpkms, 'nif'+part['name'])
			if nif_key != None:
				fpkm = fpkms[nif_key]
				part['opts']['color'] = cmapper(norm(fpkm))
			part['opts']['label'] = ''
		if part['type'] == 'Promoter':
			part['opts']['label'] = ''
		if part['type'] == 'Terminator':
			part['opts']['label'] = ''

###############################################################################
# Load all the profiles for each organism (specified in BED files)
###############################################################################

#ecoli_samples = ['N1','N2', 'EcoliMG1655_LBWB_1', 'EcoliMG1655_LBWB_2', 
#                 'EcoliMG1655_synnifI4_1', 'EcoliMG1655_synnifI4_2']
ecoli_samples = ['N2', 'EcoliMG1655_LBWB_2', 'EcoliMG1655_synnifI4_1']            
ecoli_file_prefix = '../nif_Ecoli_MG1655/results/'
ecoli_gff_prefix = '../nif_Ecoli_MG1655/data/gff/'
ecoli_gffs = load_plotting_gffs(ecoli_samples, ecoli_gff_prefix)
ecoli_fpkms = load_fpkms(ecoli_file_prefix+'fpkm.normed.matrix.txt')

#rizo_samples = ['N5','N13','N14','N15','N16','N17', 'Rhizobium_1', 'Rhizobium_2', 
#                 'Rhizobium_synnifI4_1', 'Rhizobium_synnifI4_2']
rizo_samples = ['N5','N13','N14','N15','N16','N17', 'Rhizobium_1', 'Rhizobium_2', 
                 'Rhizobium_synnifI4_1', 'Rhizobium_synnifI4_2']
rizo_file_prefix = '../nif_IRBG74_NH/results/'
rizo_gff_prefix = '../nif_IRBG74_NH/data/gff/'
rizo_gffs = load_plotting_gffs(rizo_samples, rizo_gff_prefix)
rizo_fpkms = load_fpkms(rizo_file_prefix+'fpkm.normed.matrix.txt')

koxy_samples = ['Koxym5a1_1','Koxym5a1_2', 'Koxym5a1_synnifI4_1', 
                'Koxym5a1_synnifI4_2']
koxy_file_prefix = '../nif_Koxytoca/results/'
koxy_gff_prefix = '../nif_Koxytoca/data/gff/'
koxy_gffs = load_plotting_gffs(koxy_samples, koxy_gff_prefix)
koxy_fpkms = load_fpkms(koxy_file_prefix+'fpkm.normed.matrix.txt')

pf5_samples = ['N3','N4','N18','N19','N20','N21']
pf5_file_prefix = '../nif_Pf-5/results/'
pf5_gff_prefix = '../nif_Pf-5/data/gff/'
pf5_gffs = load_plotting_gffs(pf5_samples, pf5_gff_prefix)
pf5_fpkms = load_fpkms(pf5_file_prefix+'fpkm.normed.matrix.txt')

#############################################################################
# FIG 9. FPKMs WT Cluster
#############################################################################

out_filename = 'fig_09_fpkm_wt_cluster.pdf'

max_fpkm = 50000.0
min_fpkm = 0.0

fpkms = [ecoli_fpkms['N2'], rizo_fpkms['N5'], pf5_fpkms['N3']]
design_gff = ecoli_gffs['N2']
design_chrom = 'pKU7017'
gffs = [ecoli_gffs['N2'], rizo_gffs['N5'], pf5_gffs['N3']]
chroms = ['pKU7017', 'pBBR7017', 'pPf7017']

dna_y_max = 22
dna_y_min = -22

cmapper = matplotlib.cm.get_cmap('PuBu')
cmapper = matplotlib.cm.get_cmap('RdYlBu_r')
cmapper = matplotlib.cm.get_cmap('Reds')
norm = matplotlib.colors.SymLogNorm(100, linscale=1.0, vmin=min_fpkm, vmax=max_fpkm)
#norm = matplotlib.colors.Normalize(vmin=min_fpkm, vmax=max_fpkm)

# Create the figure
fig = plt.figure(figsize=(4.4,1.5))
gs = gridspec.GridSpec(5, 1, height_ratios=[1,1,1,0.3,1])
for idx in range(len(gffs)):
	ax = plt.subplot(gs[idx])

	# Generate the design and format colour based on FPKM
	dna_design = gff_to_dnaplotlib(gffs[idx], chroms[idx])
	format_dna_design_fpkm(dna_design, fpkms[idx], cmapper, norm)
	dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
	start, end = dr.renderDNA(ax, dna_design, dr.SBOL_part_renderers())
	ax.set_xlim([start, end])
	ax.set_ylim([dna_y_min,dna_y_max])
	ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')

# Plot the DNA design
ax = plt.subplot(gs[4])
dna_design = gff_to_dnaplotlib(design_gff, design_chrom)
dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
start, end = dr.renderDNA(ax, dna_design, dr.SBOL_part_renderers())
ax.set_xlim([start, end])
ax.set_ylim([dna_y_min,dna_y_max])
ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Save the figure
plt.subplots_adjust(hspace=.00, wspace=.00, left=.05, right=.95, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True)
plt.close('all')


#############################################################################
# FIG 10. FPKMs 1.1 MB Cluster
#############################################################################

out_filename = 'fig_10_fpkm_1.1MB_cluster.pdf'

max_fpkm = 50000.0
min_fpkm = 0.0

fpkms = [ecoli_fpkms['EcoliMG1655_LBWB_2'], rizo_fpkms['N16'], pf5_fpkms['N20']]
design_gff = rizo_gffs['N16']
design_chrom = 'pBBG_MB'
gffs = [ecoli_gffs['EcoliMG1655_LBWB_2'], rizo_gffs['N16'], pf5_gffs['N20']]
chroms = ['LBWB', 'pBBG_MB', 'pPf_MB']

dna_y_max = 22
dna_y_min = -22

cmapper = matplotlib.cm.get_cmap('Reds')
norm = matplotlib.colors.SymLogNorm(100, linscale=1.0, vmin=min_fpkm, vmax=max_fpkm)
#norm = matplotlib.colors.Normalize(vmin=min_fpkm, vmax=max_fpkm)

# Create the figure
fig = plt.figure(figsize=(4.4,1.5))
gs = gridspec.GridSpec(5, 1, height_ratios=[1,1,1,0.3,1])
for idx in range(len(gffs)):
	ax = plt.subplot(gs[idx])

	# Generate the design and format colour based on FPKM
	dna_design = gff_to_dnaplotlib(gffs[idx], chroms[idx])
	format_dna_design_fpkm(dna_design, fpkms[idx], cmapper, norm)
	dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
	start, end = dr.renderDNA(ax, dna_design, dr.SBOL_part_renderers())
	ax.set_xlim([start, end])
	ax.set_ylim([dna_y_min,dna_y_max])
	ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')

# Plot the DNA design
ax = plt.subplot(gs[4])
dna_design = gff_to_dnaplotlib(design_gff, design_chrom)
dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
start, end = dr.renderDNA(ax, dna_design, dr.SBOL_part_renderers())
ax.set_xlim([start, end])
ax.set_ylim([dna_y_min,dna_y_max])
ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Save the figure
plt.subplots_adjust(hspace=.00, wspace=.00, left=.05, right=.95, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True)
plt.close('all')

#############################################################################
# FIG 11. FPKMs 2.1 i4/nif60 Cluster
#############################################################################

out_filename = 'fig_11_fpkm_2.1nif60_cluster.pdf'

max_fpkm = 50000.0
min_fpkm = 0.0

fpkms = [ecoli_fpkms['EcoliMG1655_synnifI4_1'], rizo_fpkms['N17'], pf5_fpkms['N21']]
design_gff = rizo_gffs['N17']
design_chrom = 'pBBRnif60percent_with_par'
gffs = [ecoli_gffs['EcoliMG1655_synnifI4_1'], rizo_gffs['N17'], pf5_gffs['N21']]
chroms = ['synI4', 'pBBRnif60percent_with_par', 'pPf_nif60']

dna_y_max = 22
dna_y_min = -22

cmapper = matplotlib.cm.get_cmap('Reds')
norm = matplotlib.colors.SymLogNorm(100, linscale=1.0, vmin=min_fpkm, vmax=max_fpkm)
#norm = matplotlib.colors.Normalize(vmin=min_fpkm, vmax=max_fpkm)

# Create the figure
fig = plt.figure(figsize=(4.4,1.5))
gs = gridspec.GridSpec(5, 1, height_ratios=[1,1,1,0.3,1])
for idx in range(len(gffs)):
	ax = plt.subplot(gs[idx])

	# Generate the design and format colour based on FPKM
	dna_design = gff_to_dnaplotlib(gffs[idx], chroms[idx])
	format_dna_design_fpkm(dna_design, fpkms[idx], cmapper, norm)
	dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
	start, end = dr.renderDNA(ax, dna_design, dr.SBOL_part_renderers())
	ax.set_xlim([start, end])
	ax.set_ylim([dna_y_min,dna_y_max])
	ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')

# Plot the DNA design
ax = plt.subplot(gs[4])
dna_design = gff_to_dnaplotlib(design_gff, design_chrom)
dr = dpl.DNARenderer(scale=1.0, linewidth=0.8)
start, end = dr.renderDNA(ax, dna_design, dr.SBOL_part_renderers())
ax.set_xlim([start, end])
ax.set_ylim([dna_y_min,dna_y_max])
ax.plot([start-(start*1.2),end+(end*1.2)], [0,0], color=(0,0,0), linewidth=1.2, zorder=1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Save the figure
plt.subplots_adjust(hspace=.00, wspace=.00, left=.05, right=.95, top=.99, bottom=.01)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True)
plt.close('all')








wt_gff = ecoli_gffs['N2']
wt_chrom = 'pKU7017'
v1_0_gff = rizo_gffs['N13']
v1_0_chrom = 'pBBR_OR'
v1_1_gff = rizo_gffs['N14']
v1_1_chrom = 'pBBG_MB'
v2_1_gff = rizo_gffs['N17']
v2_1_chrom = 'pBBRnif60percent_with_par'
