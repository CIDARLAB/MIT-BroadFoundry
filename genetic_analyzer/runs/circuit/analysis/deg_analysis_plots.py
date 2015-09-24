#!/usr/bin/env python
"""
DEG analysis
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
import csv

DATA_PREFIX = './fpkm_matrix/'
DEG_PREFIX = '../results/'
OUT_PREFIX = './deg_analysis/'

matplotlib.rcParams['lines.dash_joinstyle']  = 'miter'
matplotlib.rcParams['lines.dash_capstyle']   = 'butt'
matplotlib.rcParams['lines.solid_joinstyle'] = 'miter'
matplotlib.rcParams['lines.solid_capstyle']  = 'projecting'
# Make text editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype']          = 42 

all_states = [0,1,2,3,4,5,6,7]

ind_states = {}
ind_states['IPTG'] = [1,3,5,7]
ind_states['aTc']  = [2,3,6,7]
ind_states['Ara']  = [4,5,6,7]

no_ind_states = {}
no_ind_states['IPTG'] = list(set(all_states)-set(ind_states['IPTG']))
no_ind_states['aTc']  = list(set(all_states)-set(ind_states['aTc']))
no_ind_states['Ara']  = list(set(all_states)-set(ind_states['Ara']))

plt.rcParams['ytick.major.pad']='1' 
plt.rcParams['xtick.major.pad']='3'

fmt_axis_outline_width = 1.0
fmt_label_size = 8.0

###############################################################################
# PLOTTING FUNCTIONS
###############################################################################

def load_gene_fpkms (filename_in, data_cols, inc_all=False):
	gene_data = []
	circuit_data = []
	for idx in range(len(data_cols)):
		gene_data.append({})
		circuit_data.append({})

	with open(filename_in, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		# Skip header
		next(reader)
		for row in reader:
			# Classify gene
			if inc_all == False:
				if row[0].split('_')[0] == 'SYNTHETIC':
					for idx in range(len(data_cols)):
						circuit_data[idx][row[1]] = float(row[data_cols[idx]])
				else:
					for idx in range(len(data_cols)):
						gene_data[idx][row[1]] = float(row[data_cols[idx]])
			else:
				for idx in range(len(data_cols)):
						gene_data[idx][row[1]] = float(row[data_cols[idx]])
		return gene_data, circuit_data

def plot_gene_scatter (gene_data, circuit_data, unind_state, ind_states, filename_out):
	fig = plt.figure(figsize=(6.7,1.7))
	for ax_idx in range(4):
		ax = plt.subplot(1,4,ax_idx+1)
		# Gene data
		x_data = gene_data[unind_state].values()
		y_data = gene_data[ind_states[ax_idx]].values()
		ax.scatter(x_data, y_data, zorder=10, s=5, color=(0,0,1,0.1), edgecolor=(0,0,0), linewidth=0.0)
		# Circuit data
		x_data = circuit_data[unind_state].values()
		y_data = circuit_data[ind_states[ax_idx]].values()
		ax.scatter(x_data, y_data, zorder=20, s=5, color=(1,0,0,1), edgecolor=(0,0,0), linewidth=0.0)
		
		labs = circuit_data[unind_state].keys()
		for lab_idx in range(len(x_data)):
			ax.annotate(labs[lab_idx], (x_data[lab_idx],y_data[lab_idx]), fontsize=6.0, zorder=110)

		x_y_line = [0, 10000000]
		ax.plot(x_y_line, x_y_line, color=(0,0,0), linestyle='-', linewidth=0.5, zorder=100)

		line_thres_low = [1, 100000]
		line_thres_high = [10, 1000000]
		ax.plot(line_thres_low, line_thres_high, color=(0,0,0), linestyle='--', linewidth=0.5, zorder=100)
		ax.plot(line_thres_high, line_thres_low, color=(0,0,0), linestyle='--', linewidth=0.5, zorder=100)

		# Plot formatting
		ax.set_ylim([10,1000000])
		ax.set_yscale('log', linthreshy=1)
		ax.set_xlim([10,1000000])
		ax.set_xscale('log', linthreshy=1)
		ax.tick_params(axis='x', labelsize=fmt_label_size)
		ax.tick_params(axis='y', labelsize=fmt_label_size)
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
			ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.15, wspace=0.3)
	fig.savefig(filename_out, transparent=True)

def load_de_genes (filename_in, p_val=0.01, fdr=0.05):
	de_genes = {}
	with open(filename_in, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		# Skip header
		next(reader)
		for row in reader:
			if len(row) == 5:
				gene_parts = row[0].split('_locusTag_')
				if len(gene_parts) == 2:
					cur_gene = (row[0].split('_locusTag_'))[0]
					cur_p_val = float(row[3])
					cur_fdr = float(row[4])
					cur_fc = float(row[1])
					if cur_p_val <= p_val and cur_fdr <= fdr:
						de_genes[cur_gene] = [cur_p_val, cur_fdr, cur_fc]
	return de_genes

def avg_sd_exp (gene_data, states):
	new_gene_data = {}
	for k in gene_data[0].keys():
		exp_vals = []
		for s in states:
			exp_vals.append(gene_data[s][k])
		new_gene_data[k] = [np.mean(exp_vals), np.std(exp_vals)]
	return new_gene_data

def plot_inducer_scatter (gene_data, circuit_data, unind_states, ind_states, de_genes, filename_out):
	fig = plt.figure(figsize=(1.5,1.5))
	ax = plt.subplot(1,1,1)

	ind_exp = avg_sd_exp(gene_data, ind_states)
	unind_exp = avg_sd_exp(gene_data, unind_states)

	x_data = []
	y_data = []
	x_data_de = []
	y_data_de = []
	x_data_de_sd = []
	y_data_de_sd = []
	for k in ind_exp.keys():
		if k in de_genes.keys():
			x_data_de.append(unind_exp[k][0])
			x_data_de_sd.append(unind_exp[k][1])
			y_data_de.append(ind_exp[k][0])
			y_data_de_sd.append(ind_exp[k][1])
		else:
			x_data.append(unind_exp[k][0])
			y_data.append(ind_exp[k][0])

	ax.scatter(x_data, y_data, zorder=-1, s=1.1, color=(0.5,0.5,0.5), edgecolor=(0,0,0), linewidth=0.0)
	ax.scatter(x_data_de, y_data_de, zorder=20, s=5, color=(1,0,0), edgecolor=(0,0,0), linewidth=0.0)
	
	# annotate DE genes
	labs = de_genes.keys()
	for lab_idx in range(len(x_data_de)):
		ax.annotate(labs[lab_idx], (x_data_de[lab_idx],y_data_de[lab_idx]), fontsize=6.0, zorder=110)

	x_y_line = [0, 1000000]
	ax.plot(x_y_line, x_y_line, color=(0,0,0), linestyle='-', linewidth=0.5, zorder=100)

	line_thres_low = [0.1, 100000]
	line_thres_high = [1, 1000000]
	ax.plot(line_thres_low, line_thres_high, color=(0,0,0), linestyle='--', linewidth=0.5, zorder=100)
	ax.plot(line_thres_high, line_thres_low, color=(0,0,0), linestyle='--', linewidth=0.5, zorder=100)

	# Plot formatting
	ax.set_ylim([0.1,20000])
	ax.set_yscale('log', linthreshy=1)
	ax.set_xlim([0.1,20000])
	ax.set_xscale('log', linthreshy=1)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.tick_params('both', length=2, width=0.3, which='major')
	ax.tick_params('both', length=1.5, width=0.2, which='minor')
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.15, right=0.92, top=0.92, bottom=0.15, wspace=0.3)
	fig.savefig(filename_out, transparent=True)

	# Calculate R^2
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_data+x_data_de, y_data+y_data_de)
	return np.power(r_value, 2.0)

# Working
fdr_to_use = 0.01
gene_data, circuit_data = load_gene_fpkms(DATA_PREFIX + 'fpkm_data.txt', [2,3,4,5,6,7,8,9])
de_genes_Ara = load_de_genes (DEG_PREFIX+'ara_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_IPTG = load_de_genes (DEG_PREFIX+'iptg_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_aTc = load_de_genes (DEG_PREFIX+'atc_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)

Ara_R_sq = plot_inducer_scatter(gene_data, circuit_data, no_ind_states['Ara'], ind_states['Ara'], de_genes_Ara, OUT_PREFIX+'scatter_Ara.pdf')
IPTG_R_sq = plot_inducer_scatter(gene_data, circuit_data, no_ind_states['IPTG'], ind_states['IPTG'], de_genes_IPTG, OUT_PREFIX+'scatter_IPTG.pdf')
aTc_R_sq = plot_inducer_scatter(gene_data, circuit_data, no_ind_states['aTc'], ind_states['aTc'], de_genes_aTc, OUT_PREFIX+'scatter_aTc.pdf')
print('R-squared: Ara='+str(Ara_R_sq)+', IPTG='+str(IPTG_R_sq)+', aTc='+str(aTc_R_sq))

# Compare the all and broken (6,7,8) states
de_genes_flask_vs_tube_all = load_de_genes (DEG_PREFIX+'flask_vs_tube.de.analysis.txt', p_val=0.01, fdr=1.0)
de_genes_flask_vs_tube_678 = load_de_genes (DEG_PREFIX+'broken_flask_vs_tube.de.analysis.txt', p_val=0.01, fdr=1.0)

def split_de_gene_up_down (de_genes):
	up_genes = []
	down_genes = []
	for g in de_genes.keys():
		if de_genes[g][2] < 0:
			down_genes.append(g)
		else:
			up_genes.append(g)
	return set(up_genes), set(down_genes)

de_genes_flask_vs_tube_all_up, de_genes_flask_vs_tube_all_down = split_de_gene_up_down(de_genes_flask_vs_tube_all)
de_genes_flask_vs_tube_678_up, de_genes_flask_vs_tube_678_down = split_de_gene_up_down(de_genes_flask_vs_tube_678)

print('UP REGULATED:')
print('Only All: '+str(len(de_genes_flask_vs_tube_all_up.difference(de_genes_flask_vs_tube_678_up))))
print('Shared: '+str(len(de_genes_flask_vs_tube_all_up.intersection(de_genes_flask_vs_tube_678_up))))
print('Only 6,7,8: '+str(len(de_genes_flask_vs_tube_678_up.difference(de_genes_flask_vs_tube_all_up))))

print('\nDOWN REGULATED:')
print('Only All: '+str(len(de_genes_flask_vs_tube_all_down.difference(de_genes_flask_vs_tube_678_down))))
print('Shared: '+str(len(de_genes_flask_vs_tube_all_down.intersection(de_genes_flask_vs_tube_678_down))))
print('Only 6,7,8: '+str(len(de_genes_flask_vs_tube_678_down.difference(de_genes_flask_vs_tube_all_down))))


########################################################
# FOLD-CHANGE MATRIX FOR TRNASPORTERS
########################################################

def fc_matrix (fpkm_data_1, fpkm_data_2, genes):
	fc_data = np.zeros( (len(fpkm_data_1), len(genes)) )
	for s in range(len(fpkm_data_1)):
		for g_idx in range(len(genes)):
			g = genes[g_idx]
			data_1 = fpkm_data_1[s][g]
			data_2 = fpkm_data_2[s][g]
			fc = 0.0
			if data_1 != 0.0 :
				fc = np.log2(data_2/data_1)
			else:
				fc = 100.0
			fc_data[s,g_idx] = fc
	return np.array(fc_data)

genes_to_inc = ['AraC', 'araE', 'araF', 'araG', 'araH', 'malE', 'malF', 'malG', 'malK', 'xylA', 'xylB', 'xylF', 'xylG', 'xylH', 'xylR']
gene_data_1, circuit_data_1 = load_gene_fpkms(DATA_PREFIX + 'fpkm_data.txt', [2,3,4,5,6,7,8,9], inc_all=True)
gene_data_2, circuit_data_2 = load_gene_fpkms(DATA_PREFIX + 'fpkm_data.txt', [10,11,12,13,14,15,16,17], inc_all=True)
fc_matrix_data = fc_matrix (gene_data_1, gene_data_2, genes_to_inc)

input_states = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],
                [0,0,1],[1,0,1],[0,1,1],[1,1,1]]

fig = plt.figure(figsize=(3.5,1.3))
ax_input = fig.add_axes([0.02,0.02,0.13,0.96])
ax_input.matshow(input_states, aspect='auto', cmap=plt.cm.Greys)
ax = fig.add_axes([0.15,0.02,0.7,0.96])
im = ax.matshow(fc_matrix_data, aspect='auto', origin='upper', cmap=plt.cm.RdBu_r, vmin=-4.0, vmax=4.0)

pcolor = fig.add_axes([0.88,0.02,0.015,0.96])
cbar = plt.colorbar(im, cax=pcolor)
cbar.locator = matplotlib.ticker.FixedLocator([])
cbar.update_ticks()

for axis in ['top','right', 'left', 'bottom']:
	ax_input.spines[axis].set_linewidth(0.8)
	ax.spines[axis].set_linewidth(0.8)
cbar.ax.get_children()[2].set_linewidth(0.8)
ax_input.set_yticks([])
ax_input.set_xticks([])
ax.set_yticks([])
ax.set_xticks([])

fig.savefig(OUT_PREFIX+'transporter_fc_heatmap.pdf', transparent=True)
plt.close('all')
