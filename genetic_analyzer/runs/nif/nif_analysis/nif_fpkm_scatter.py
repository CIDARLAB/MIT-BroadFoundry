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

plt.rcParams['ytick.major.pad']='1'
plt.rcParams['xtick.major.pad']='3'
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


def format_scatter_axes (ax):
	fmt_label_size = 11.0
	ax.tick_params(axis='y', which='major', labelsize=fmt_label_size, pad=1, length=2.5, width=0.5)
	ax.tick_params(axis='x', which='major', labelsize=fmt_label_size, pad=3, length=2.5, width=0.5)
	ax.set_xscale('symlog',linthreshy=1.1)
	ax.set_yscale('symlog',linthreshy=1.1)
	ax.set_xlim([0,300000])
	ax.set_ylim([0,300000])
	ax.plot([0,10000000],[0,10000000], linestyle='--', linewidth=1.0, color=(0.5,0.5,0.5), zorder=-10)

def plot_fpkm_scatter_wt (ax, fpkms_wt, fpkms_y, genes, wt_genes):
	# extract the data
	x_vals = {}
	for tag in fpkms_wt.keys():
		if tag in wt_genes:
			use_tag = tag.split('_')[0]
			#print tag, use_tag
			x_vals[use_tag] = fpkms_wt[tag]
	y_vals = {}
	for tag in fpkms_y.keys():
		if tag in genes:
			y_vals[tag] = fpkms_y[tag]
	x_points = []
	y_points = []
	for g in genes:
		x_points.append(x_vals[g])
		y_points.append(y_vals[g])
	# plot the scatter
	ax.scatter(x_points, y_points, s=30, facecolor=(1,1,1,0), linewidth=1.2)
	format_scatter_axes(ax)

def plot_fpkm_scatter (ax, fpkms_x, fpkms_y, genes):
	# extract the data
	x_vals = {}
	for tag in fpkms_x.keys():
		if tag in genes:
			x_vals[tag] = fpkms_x[tag]
	y_vals = {}
	for tag in fpkms_y.keys():
		if tag in genes:
			y_vals[tag] = fpkms_y[tag]
	x_points = []
	y_points = []
	for g in genes:
		x_points.append(x_vals[g])
		y_points.append(y_vals[g])
	# plot the scatter
	ax.scatter(x_points, y_points, s=30, facecolor=(1,1,1,0), linewidth=1.5)
	format_scatter_axes(ax)

###############################################################################
# Load all the profiles for each organism (specified in BED files)
###############################################################################

ecoli_samples = ['N1','N2', 'EcoliMG1655_LBWB_1', 'EcoliMG1655_LBWB_2', 
                 'EcoliMG1655_synnifI4_1', 'EcoliMG1655_synnifI4_2']         
ecoli_file_prefix = '../nif_Ecoli_MG1655/results/'
ecoli_fpkms = load_fpkms(ecoli_file_prefix+'fpkm.normed.matrix.txt')

rizo_samples = ['N5','N13','N14','N15','N16','N17', 'Rhizobium_1', 'Rhizobium_2', 
                 'Rhizobium_synnifI4_1', 'Rhizobium_synnifI4_2']
rizo_file_prefix = '../nif_IRBG74_NH/results/'
rizo_fpkms = load_fpkms(rizo_file_prefix+'fpkm.normed.matrix.txt')

koxy_samples = ['Koxym5a1_1','Koxym5a1_2', 'Koxym5a1_synnifI4_1', 
                'Koxym5a1_synnifI4_2']
koxy_file_prefix = '../nif_Koxytoca/results/'
koxy_fpkms = load_fpkms(koxy_file_prefix+'fpkm.normed.matrix.txt')

pf5_samples = ['N3','N4','N18','N19','N20','N21']
pf5_file_prefix = '../nif_Pf-5/results/'
pf5_fpkms = load_fpkms(pf5_file_prefix+'fpkm.normed.matrix.txt')

#############################################################################
# FIG 9. Scatter plots
#############################################################################

out_filename = 'fig_09_fpkm_wt_cluster.1MB_cluster.pdf'

max_fpkm = 50000.0
min_fpkm = 0.0

genes_wt_kleb = ['nifA_KOXM_26810', 'nifB_DBG_0002', 'nifD_KOXM_26890', 'nifE_DBG_0001', 
                 'nifF_KOXM_26820', 'nifH_KOXM_26895', 'nifJ_KOXM_26900', 'nifK_KOXM_26885', 
                 'nifL_KOXM_26815', 'nifM_KOXM_26825', 'nifN_KOXM_26860', 'nifQ_KOXM_26795', 
                 'nifS_KOXM_26845', 'nifT_KOXM_26880', 'nifU_KOXM_26850', 'nifV_KOXM_26840', 
                 'nifW_KOXM_26835', 'nifX_KOXM_26855', 'nifY_KOXM_26875', 'nifZ_KOXM_26830']

genes_wt = ['nifU', 'nifS', 'nifV', 'nifW', 'nifH', 'nifD', 'nifK', 'nifT', 
            'nifY', 'nifE', 'nifN', 'nifX', 'nifA', 'nifL', 'nifQ', 'nifB',
            'nifJ', 'nifF']

genes_refactored = ['nifH', 'nifD', 'nifK', 'nifY', 'nifE', 'nifN', 'nifJ', 
                    'nifB', 'nifQ', 'nifF', 'nifU', 'nifS', 'nifV', 'nifW', 
                    'nifZ', 'nifM',]

fpkms_wt = [koxy_fpkms['Koxym5a1_1'], 
            ecoli_fpkms['N2'], 
            rizo_fpkms['N5'],
            pf5_fpkms['N3']]

fpkms_v1_1 = [koxy_fpkms['Koxym5a1_synnifI4_2'], 
              ecoli_fpkms['EcoliMG1655_LBWB_1'], 
              rizo_fpkms['N16'],
              pf5_fpkms['N20']]

fpkms_v2_1 = [koxy_fpkms['Koxym5a1_synnifI4_2'], 
              ecoli_fpkms['EcoliMG1655_synnifI4_1'], 
              rizo_fpkms['N17'],
              pf5_fpkms['N21']]

# Create the figure
fig = plt.figure(figsize=(7.0,4.5))
gs = gridspec.GridSpec(2, 3, width_ratios=[1,1,1], height_ratios=[1,1])

ax = plt.subplot(gs[0])
plot_fpkm_scatter_wt (ax, fpkms_wt[0], fpkms_wt[1], genes_wt, genes_wt_kleb)
ax = plt.subplot(gs[1])
plot_fpkm_scatter_wt (ax, fpkms_wt[0], fpkms_wt[2], genes_wt, genes_wt_kleb)
ax = plt.subplot(gs[2])
plot_fpkm_scatter_wt (ax, fpkms_wt[0], fpkms_wt[3], genes_wt, genes_wt_kleb)

ax = plt.subplot(gs[3])
plot_fpkm_scatter (ax, fpkms_v2_1[0], fpkms_v2_1[1], genes_refactored)
ax = plt.subplot(gs[4])
plot_fpkm_scatter (ax, fpkms_v2_1[0], fpkms_v2_1[2], genes_refactored)
ax = plt.subplot(gs[5])
plot_fpkm_scatter (ax, fpkms_v2_1[0], fpkms_v2_1[3], genes_refactored)

#ax = plt.subplot(gs[6])
#plot_fpkm_scatter (ax, fpkms_v2_1[0], fpkms_v2_1[1], genes_refactored)
#ax = plt.subplot(gs[7])
#plot_fpkm_scatter (ax, fpkms_v2_1[0], fpkms_v2_1[2], genes_refactored)
#ax = plt.subplot(gs[8])
#plot_fpkm_scatter (ax, fpkms_v2_1[0], fpkms_v2_1[3], genes_refactored)

# Save the figure
plt.subplots_adjust(hspace=.20, wspace=.20, left=.05, right=.99, top=.99, bottom=.05)
fig.savefig(OUTPUT_PREFIX+out_filename, transparent=True, dpi=600)
plt.close('all')

