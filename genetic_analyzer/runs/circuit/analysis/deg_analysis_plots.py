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

on_states = {}
on_states['AmtR']  = [1,2,3,5,6,7]
on_states['LitR']  = [2,3,4,5,6,7]
on_states['BM3R1'] = [4,5,6,7]
on_states['SrpR']  = [0,1,2,3,4]
on_states['PhlF']  = [0,1,5,6,7]

off_states = {}
off_states['AmtR'] = list(set(all_states)-set(on_states['AmtR']))
off_states['LitR'] = list(set(all_states)-set(on_states['LitR']))
off_states['BM3R1'] = list(set(all_states)-set(on_states['BM3R1']))
off_states['SrpR'] = list(set(all_states)-set(on_states['SrpR']))
off_states['PhlF'] = list(set(all_states)-set(on_states['PhlF']))

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
	#fig = plt.figure(figsize=(1.0,1.0))
	fig = plt.figure(figsize=(1.1,1.1))
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

	# 1.1, 3
	ax.scatter(x_data, y_data, zorder=-1, s=0.6, color=(0.5,0.5,0.5), edgecolor=(0,0,0), linewidth=0.0)
	ax.scatter(x_data_de, y_data_de, zorder=20, s=1.2, color=(1,0,0), edgecolor=(0,0,0), linewidth=0.0)
	
	# annotate DE genes
	#labs = de_genes.keys()
	#for lab_idx in range(len(x_data_de)):
	#	ax.annotate(labs[lab_idx], (x_data_de[lab_idx],y_data_de[lab_idx]), fontsize=6.0, zorder=110)

	x_y_line = [0, 1000000]
	ax.plot(x_y_line, x_y_line, color=(0,0,0), linestyle='-', linewidth=0.6, zorder=100)

	#line_thres_low = [0.1, 100000]
	#line_thres_high = [1, 1000000]
	#ax.plot(line_thres_low, line_thres_high, color=(0,0,0), linestyle='--', linewidth=0.5, zorder=100)
	#ax.plot(line_thres_high, line_thres_low, color=(0,0,0), linestyle='--', linewidth=0.5, zorder=100)

	# Plot formatting
	ax.set_ylim([0.5,30000])
	ax.set_yscale('log', linthreshy=1)
	ax.set_xlim([0.5,30000])
	ax.set_xscale('log', linthreshy=1)
	ax.tick_params(axis='x', labelsize=7.0)
	ax.tick_params(axis='y', labelsize=7.0)
	ax.tick_params('both', length=2, width=0.3, which='major')
	ax.tick_params('both', length=1.2, width=0.15, which='minor')

	ax.tick_params(axis='y', which='major', pad=0)
	ax.tick_params(axis='x', which='major', pad=2)

	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(0.7)
		ax.spines[axis].set_linewidth(0.7)
		ax.spines[axis].set_linewidth(0.7)
	plt.subplots_adjust(left=0.18, right=0.99, top=0.99, bottom=0.15)
	fig.savefig(filename_out, transparent=True)

	# Calculate R^2
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_data+x_data_de, y_data+y_data_de)
	return np.power(r_value, 2.0)

# Working
fdr_to_use = 0.05
gene_data, circuit_data = load_gene_fpkms(DATA_PREFIX + 'fpkm_data.txt', [2,3,4,5,6,7,8,9])
de_genes_Ara = load_de_genes (DEG_PREFIX+'ara_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_IPTG = load_de_genes (DEG_PREFIX+'iptg_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_aTc = load_de_genes (DEG_PREFIX+'atc_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)

de_genes_AmtR = load_de_genes (DEG_PREFIX+'AmtR_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_LitR = load_de_genes (DEG_PREFIX+'LitR_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_BM3R1 = load_de_genes (DEG_PREFIX+'BM3R1_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_SrpR = load_de_genes (DEG_PREFIX+'SrpR_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
de_genes_PhlF = load_de_genes (DEG_PREFIX+'PhlF_comp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)

print de_genes_PhlF

Ara_R_sq = plot_inducer_scatter(gene_data, circuit_data, no_ind_states['Ara'], ind_states['Ara'], de_genes_Ara, OUT_PREFIX+'scatter_Ara.pdf')
IPTG_R_sq = plot_inducer_scatter(gene_data, circuit_data, no_ind_states['IPTG'], ind_states['IPTG'], de_genes_IPTG, OUT_PREFIX+'scatter_IPTG.pdf')
aTc_R_sq = plot_inducer_scatter(gene_data, circuit_data, no_ind_states['aTc'], ind_states['aTc'], de_genes_aTc, OUT_PREFIX+'scatter_aTc.pdf')
print('R-squared: Ara='+str(Ara_R_sq)+', IPTG='+str(IPTG_R_sq)+', aTc='+str(aTc_R_sq))

AmtR_R_sq = plot_inducer_scatter(gene_data, circuit_data, off_states['AmtR'] , on_states['AmtR'] , de_genes_AmtR, OUT_PREFIX+'scatter_AmtR.pdf')
LitR_R_sq = plot_inducer_scatter(gene_data, circuit_data, off_states['LitR'] , on_states['LitR'] , de_genes_LitR, OUT_PREFIX+'scatter_LitR.pdf')
BM3R1_R_sq = plot_inducer_scatter(gene_data, circuit_data, off_states['BM3R1'] , on_states['BM3R1'] , de_genes_BM3R1, OUT_PREFIX+'scatter_BM3R1.pdf')
SrpR_R_sq = plot_inducer_scatter(gene_data, circuit_data, off_states['SrpR'] , on_states['SrpR'] , de_genes_SrpR, OUT_PREFIX+'scatter_SrpR.pdf')
PhlF_R_sq = plot_inducer_scatter(gene_data, circuit_data, off_states['PhlF'] , on_states['PhlF'] , de_genes_PhlF, OUT_PREFIX+'scatter_PhlF.pdf')
print('R-squared: AmtR='+str(AmtR_R_sq)+', LitR='+str(LitR_R_sq)+', BM3R1='+str(BM3R1_R_sq)+', SrpR='+str(SrpR_R_sq)+', PhlF='+str(PhlF_R_sq))


fdr_to_use = 1.0
rep_2_sq = plot_inducer_scatter(gene_data, circuit_data, [0] , [1], {}, OUT_PREFIX+'scatter_2_reps.pdf')
print('R-squared: 2 repressors vs 3 repressors='+str(rep_2_sq))
de_genes_rep_3_yfp = load_de_genes (DEG_PREFIX+'rep_3_yfp_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
rep_3_yfp_sq = plot_inducer_scatter(gene_data, circuit_data, [0,1] , [2,3,4], de_genes_rep_3_yfp, OUT_PREFIX+'scatter_3_yfp_reps.pdf')
print('R-squared: 3 repressors + YFP='+str(rep_3_yfp_sq))
de_genes_rep_4 = load_de_genes (DEG_PREFIX+'rep_4_tube.de.analysis.txt', p_val=0.01, fdr=fdr_to_use)
rep_4_sq = plot_inducer_scatter(gene_data, circuit_data, [0,1] , [5,6,7], de_genes_rep_4, OUT_PREFIX+'scatter_4_reps.pdf')
print('R-squared: 4 repressors='+str(rep_4_sq))



de_genes_State_1 = load_de_genes (DEG_PREFIX+'state_1_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_1_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [0], de_genes_State_1, OUT_PREFIX+'scatter_state_1.pdf')
print('R-squared: State 1='+str(State_1_sq))
de_genes_State_2 = load_de_genes (DEG_PREFIX+'state_2_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_2_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [1], de_genes_State_2, OUT_PREFIX+'scatter_state_2.pdf')
print('R-squared: State 2='+str(State_2_sq))
de_genes_State_3 = load_de_genes (DEG_PREFIX+'state_3_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_3_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [2], de_genes_State_3, OUT_PREFIX+'scatter_state_3.pdf')
print('R-squared: State 3='+str(State_3_sq))
de_genes_State_4 = load_de_genes (DEG_PREFIX+'state_4_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_4_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [3], de_genes_State_4, OUT_PREFIX+'scatter_state_4.pdf')
print('R-squared: State 4='+str(State_4_sq))
de_genes_State_5 = load_de_genes (DEG_PREFIX+'state_5_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_5_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [4], de_genes_State_5, OUT_PREFIX+'scatter_state_5.pdf')
print('R-squared: State 5='+str(State_5_sq))
de_genes_State_6 = load_de_genes (DEG_PREFIX+'state_6_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_6_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [5], de_genes_State_6, OUT_PREFIX+'scatter_state_6.pdf')
print('R-squared: State 6='+str(State_6_sq))
de_genes_State_7 = load_de_genes (DEG_PREFIX+'state_7_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_7_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [6], de_genes_State_7, OUT_PREFIX+'scatter_state_7.pdf')
print('R-squared: State 7='+str(State_7_sq))
de_genes_State_8 = load_de_genes (DEG_PREFIX+'state_8_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
State_8_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [7], de_genes_State_8, OUT_PREFIX+'scatter_state_8.pdf')
print('R-squared: State 8='+str(State_8_sq))

"""
de_genes_rep_2 = load_de_genes (DEG_PREFIX+'rep_2_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
rep_2_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [0], de_genes_rep_2, OUT_PREFIX+'scatter_2_reps.pdf')
print('R-squared: 2 repressors='+str(rep_2_sq))
de_genes_rep_3 = load_de_genes (DEG_PREFIX+'rep_2_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
rep_3_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [1], de_genes_rep_3, OUT_PREFIX+'scatter_3_reps.pdf')
print('R-squared: 3 repressors='+str(rep_3_sq))
de_genes_rep_3_yfp = load_de_genes (DEG_PREFIX+'rep_2_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
rep_3_yfp_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [2,3,4], de_genes_rep_3_yfp, OUT_PREFIX+'scatter_3_yfp_reps.pdf')
print('R-squared: 3 repressors + YFP='+str(rep_3_yfp_sq))
de_genes_rep_4 = load_de_genes (DEG_PREFIX+'rep_2_tube.de.analysis.txt', p_val=0.05, fdr=fdr_to_use)
rep_4_sq = plot_inducer_scatter(gene_data, circuit_data, all_states , [5,6,7], de_genes_rep_4, OUT_PREFIX+'scatter_4_reps.pdf')
print('R-squared: 4 repressors='+str(rep_4_sq))
"""



def print_de_gene_stats (name, de_genes):
	num_down = 0
	num_up = 0
	for k in de_genes.keys():
		g = de_genes[k]
		if g[2] > 0.0:
			num_down += 1
		else:
			num_up += 1
	print name, 'UP:', num_up, 'DOWN:', num_down

print_de_gene_stats('Ara', de_genes_Ara)
print_de_gene_stats('aTc', de_genes_aTc)

print_de_gene_stats('AmtR', de_genes_AmtR)
print_de_gene_stats('LitR', de_genes_LitR)
print_de_gene_stats('BM3R1', de_genes_BM3R1)
print_de_gene_stats('SrpR', de_genes_SrpR)
print_de_gene_stats('PhlF', de_genes_PhlF)

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

def make_study (de_genes, filename_out_up, filename_out_down):
	f_out_up = open(filename_out_up, 'w')
	f_out_down = open(filename_out_down, 'w')
	for k in de_genes.keys():
		g = de_genes[k]
		if g[2] > 0.0:
			f_out_down.write(k+'\n')
		else:
			f_out_up.write(k+'\n')
	f_out_up.close()
	f_out_down.close()

STUDY_PREFIX = './deg_analysis/studies/'
make_study(de_genes_Ara, STUDY_PREFIX+'study_tube_Ara_up.txt', STUDY_PREFIX+'study_tube_Ara_down.txt')
make_study(de_genes_aTc, STUDY_PREFIX+'study_tube_aTc_up.txt', STUDY_PREFIX+'study_tube_aTc_down.txt')
make_study(de_genes_IPTG, STUDY_PREFIX+'study_tube_IPTG_up.txt', STUDY_PREFIX+'study_tube_IPTG_down.txt')

make_study(de_genes_AmtR, STUDY_PREFIX+'study_tube_AmtR_up.txt', STUDY_PREFIX+'study_tube_AmtR_down.txt')
make_study(de_genes_LitR, STUDY_PREFIX+'study_tube_LitR_up.txt', STUDY_PREFIX+'study_tube_LitR_down.txt')
make_study(de_genes_BM3R1, STUDY_PREFIX+'study_tube_BM3R1_up.txt', STUDY_PREFIX+'study_tube_BM3R1_down.txt')
make_study(de_genes_SrpR, STUDY_PREFIX+'study_tube_SrpR_up.txt', STUDY_PREFIX+'study_tube_SrpR_down.txt')
make_study(de_genes_PhlF, STUDY_PREFIX+'study_tube_PhlF_up.txt', STUDY_PREFIX+'study_tube_PhlF_down.txt')

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
#ax_input = fig.add_axes([0.02,0.02,0.13,0.96])
#ax_input.matshow(input_states, aspect='auto', cmap=plt.cm.Greys)
ax = fig.add_axes([0.15,0.02,0.7,0.96])
im = ax.matshow(fc_matrix_data, aspect='auto', origin='upper', cmap=plt.cm.RdBu_r, vmin=-4.0, vmax=4.0)

pcolor = fig.add_axes([0.88,0.02,0.015,0.96])
cbar = plt.colorbar(im, cax=pcolor)
cbar.locator = matplotlib.ticker.FixedLocator([])
cbar.update_ticks()

for axis in ['top','right', 'left', 'bottom']:
	#ax_input.spines[axis].set_linewidth(0.8)
	ax.spines[axis].set_linewidth(0.8)
cbar.ax.get_children()[2].set_linewidth(0.8)
#ax_input.set_yticks([])
#ax_input.set_xticks([])
ax.set_yticks([])
ax.set_xticks([])

fig.savefig(OUT_PREFIX+'transporter_fc_heatmap.pdf', transparent=True)
plt.close('all')
