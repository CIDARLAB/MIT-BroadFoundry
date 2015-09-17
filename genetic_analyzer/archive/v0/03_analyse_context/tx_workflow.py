#!/usr/bin/env python
"""Full transcriptional analysis workflow
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

from workflow_utils import *
import promoter_classifier as pc
import terminator_classifier as tc

import gene_cluster_library as gcl
import gene_cluster_analysis as gca
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

INPUT_DIR = './_data/'
RESULTS_DIR = './_results/'
SEQ_RESULTS_DIR = RESULTS_DIR + '_seq_analysis/'

PLOT_PROMOTER_ARCH = True
PLOT_TERMINATOR_ARCH = True

USE_TRACE = False

###############################################################################
# LOAD DATA, CALCULATE STRENGTHS
###############################################################################

# Load the Stata nif library data
nifs = gcl.GeneClusterLibrary()
nifs.load(INPUT_DIR+'clean_nif_stata_library.txt')

# Classify all promoter and terminator contexts in the library
p_classifier = pc.PromoterClassifier()
t_classifier = tc.TerminatorClassifier()
p_data, p_bad_data = p_classifier.classify(nifs)
t_data, t_bad_data = t_classifier.classify(nifs)

# Load the RNA-seq nt-level data
tx_traces = gca.load_stata_strand_data(INPUT_DIR+'phys_depths3.csv')

# Promoter strength calculations
if USE_TRACE == True:
	p_trace_strength = pc.PromoterStrengthFromTrace()
	p_trace_strength.append_strength(nifs, tx_traces, p_data)
else:
	p_rsem_strength = pc.PromoterStrengthFromRSEM(INPUT_DIR+'promoter_express_pref.txt')
	p_rsem_strength.append_strength(p_data)

# Terminator strength calculations
if USE_TRACE == True:
	t_trace_strength = tc.TerminationFromTrace(fwd_skip=0, fwd_len=15, rev_skip=0, rev_len=15)
	t_trace_strength.append_strength(nifs, tx_traces, t_data)
else:
	t_rsem_strength = tc.TerminationFromRSEM(INPUT_DIR+'terminator_express_pref.txt')
	t_rsem_strength.append_strength(t_data)

# Remove some of the axis lines
def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

###############################################################################
# 1. OVERALL DISTRIBUTIONS
###############################################################################

def prom_PA_PP_detail_plot (p_data, key, p_contexts, prom_bin_list, x_range, y_max, p_color, y_plots, out_filename):
	fig = plt.figure(figsize=(3.5,y_plots*2))
	ax_i = 1
	p_data_contexts = split_on_detail(p_data, key)
	for c in p_contexts:
		ax = plt.subplot(len(p_contexts),1,ax_i)
		simpleaxis(ax)
		if c in p_data_contexts.keys():
			cur_data = [x for x in all_strengths(p_data_contexts[c]) if x != None]
			ax.hist(cur_data, bins=prom_bin_list, color=p_color, edgecolor=p_color)
			plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)
			if np.median(cur_data) > 0.0:
				plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
			else:
				plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
			plt.text(0.02, 0.95, plot_stats, fontsize=14,
		 		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
			
		plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
		#ax.set_title(c, fontsize=8)
		ax.tick_params(axis='x', labelsize=14)
		ax.tick_params(axis='y', labelsize=14)
		ax.set_xlim([x_range[0],x_range[1]])
		ax.set_ylim([0,y_max])
		if ax_i < len(p_contexts):
			ax.set_xticklabels([])
			ax.set_axisbelow(True)
		ax_i += 1
	plt.tight_layout()
	plt.savefig(out_filename)

def prom_PA_plot (p_data_contexts, p_contexts, prom_bin_list, p_color, x_range, y_max, out_prefix):
	fig = plt.figure(figsize=(3.5,12))
	ax_i = 1
	for c in p_contexts:
		ax = plt.subplot(len(p_contexts),1,ax_i)
		simpleaxis(ax)
		if c in p_data_contexts.keys():
			cur_data = [x for x in all_strengths(p_data_contexts[c]) if x != None]
			ax.hist(cur_data, bins=prom_bin_list, color=p_color, edgecolor=p_color)
			plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)

			if np.median(cur_data) > 0.0:
				plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
			else:
				plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
			plt.text(0.02, 0.95, plot_stats, fontsize=14,
		 		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
			
			nifCDS_names = ['nifW', 'nifZ', 'nifM', 'nifU', 'nifV', 'nifS']
			prom_PA_PP_detail_plot(p_data_contexts[c], 'ds_CDS_name', nifCDS_names, prom_bin_list, x_range, 5, p_color, len(nifCDS_names), out_prefix+'_'+c+'_ds_CDS.pdf')
			prom_PA_PP_detail_plot(p_data_contexts[c], 'us_CDS_name', nifCDS_names, prom_bin_list, x_range, 5, p_color, len(nifCDS_names), out_prefix+'_'+c+'_us_CDS.pdf')
			
			prom_names = ['P1', 'P2', 'P3']
			prom_PA_PP_detail_plot(p_data_contexts[c], 'ds_Promoter_name', prom_names, prom_bin_list, x_range, 30, p_color, len(prom_names), out_prefix+'_'+c+'_ds_Promoter.pdf')
			prom_PA_PP_detail_plot(p_data_contexts[c], 'us_Promoter_name', prom_names, prom_bin_list, x_range, 30, p_color, len(prom_names), out_prefix+'_'+c+'_us_Promoter.pdf')

			plt.figure(fig.number)
		plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
		#ax.set_title(c, fontsize=8)
		ax.tick_params(axis='x', labelsize=14)
		ax.tick_params(axis='y', labelsize=14)
		ax.set_xlim([x_range[0],x_range[1]])
		ax.set_ylim([0,y_max])
		if ax_i < len(p_contexts):
			ax.set_xticklabels([])
			ax.set_axisbelow(True)
		ax_i += 1
	plt.tight_layout()
	plt.savefig(out_prefix+'.pdf')

if PLOT_TERMINATOR_ARCH == True:
	# Terminators we do manually
	fig = plt.figure(figsize=(5,12))
	ax_i = 1
	x_min = 0.0
	x_max = 1.0
	bin_size = 0.05
	#x_min = -120
	#x_max = 120
	#bin_size = 10
	term_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
	t_data_contexts = split_on_context(t_data)
	t_contexts = ['CDS-T', 'CDS-T-P', 'CDS-T-CDS', 'CDS-T-rP', 'CDS-T-rCDS', 'CDS-T-rT']
	for c in t_contexts:
		ax = plt.subplot(len(all_contexts(p_data)),1,ax_i)
		simpleaxis(ax)
		if c in t_data_contexts.keys():
			cur_data = [x for x in all_strengths(t_data_contexts[c]) if x != None]
			ax.hist(cur_data, bins=term_bin_list, color=[0.5,0.5,0.5], edgecolor=[0.5,0.5,0.5])
			plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)
			
			if np.median(cur_data) > 0.0:
				plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
			else:
				plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
			plt.text(0.02, 0.95, plot_stats, fontsize=14,
	     		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

			nifCDS_names = ['nifW', 'nifZ', 'nifM', 'nifU', 'nifV', 'nifS']
			prom_PA_PP_detail_plot(t_data_contexts[c], 'us_CDS_name', nifCDS_names, term_bin_list, [x_min,x_max], 11, [0.5,0.5,0.5], len(nifCDS_names), RESULTS_DIR+'term_PA_contexts_'+c+'_us_CDS.pdf')
			prom_PA_PP_detail_plot(t_data_contexts[c], 'ds_CDS_name', nifCDS_names, term_bin_list, [x_min,x_max], 11, [0.5,0.5,0.5], len(nifCDS_names), RESULTS_DIR+'term_PA_contexts_'+c+'_ds_CDS.pdf')		
			prom_names = ['P1', 'P2', 'P3']
			prom_PA_PP_detail_plot(t_data_contexts[c], 'ds_Promoter_name', prom_names, term_bin_list, [x_min,x_max], 30, [0.5,0.5,0.5], len(prom_names), RESULTS_DIR+'term_PA_contexts_'+c+'_ds_Promoter.pdf')

			plt.figure(fig.number)
		plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
		ax.tick_params(axis='x', labelsize=14)
		ax.tick_params(axis='y', labelsize=14)
		ax.set_xlim([x_min,x_max])
		ax.set_ylim([0,27])
		if ax_i < 6:
			ax.set_xticklabels([])
			ax.set_axisbelow(True)
		ax_i += 1
	plt.tight_layout()
	plt.subplots_adjust(hspace=.08)
	plt.savefig(RESULTS_DIR+'term_PA_contexts.pdf')

	# All terminators plot
	fig = plt.figure(figsize=(5,3))
	ax = plt.subplot(1,1,1)
	simpleaxis(ax)
	cur_data = [x for x in all_strengths(t_data) if x != None]
	ax.hist(cur_data, bins=term_bin_list, color=[0.5,0.5,0.5], edgecolor=[0.5,0.5,0.5])
	plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)
	if np.median(cur_data) > 0.0:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	else:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	plt.text(0.02, 0.95, plot_stats, fontsize=14,
 		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
	plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.set_xlim([x_min,x_max])
	ax.set_ylim([0,50])
	plt.tight_layout()
	plt.savefig(RESULTS_DIR+'term_All.pdf')

# Parameters for the plot
x_min = 0.0
x_max = 100000
x_range = [x_min, x_max]
bin_size = 10000
#x_min = -40
#x_max = 70
#x_range = [x_min, x_max]
#bin_size = 5
prom_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
p_contexts = ['P-CDS', 'T-P-CDS', 'P-P-CDS', 'T-P-P', 'CDS-P-CDS', 'rP-P-CDS']

if PLOT_PROMOTER_ARCH == True:
	prom_PA_plot (split_on_context(filter_part_names(p_data, ['P1'])), 
		          p_contexts, prom_bin_list, [0.9,0.6,0.6], x_range, 27, 
		          RESULTS_DIR+'prom_P1_PA_contexts')

	prom_PA_plot (split_on_context(filter_part_names(p_data, ['P2'])), 
		          p_contexts, prom_bin_list, [0.6,0.9,0.6], x_range, 11, 
		          RESULTS_DIR+'prom_P2_PA_contexts')

	prom_PA_plot (split_on_context(filter_part_names(p_data, ['P3'])), 
		          p_contexts, prom_bin_list, [0.6,0.6,0.9], x_range, 11, 
		          RESULTS_DIR+'prom_P3_PA_contexts')

	# All promoters plot
	# P1
	fig = plt.figure(figsize=(5,8))
	ax = plt.subplot(3,1,1)
	simpleaxis(ax)
	cur_data = [x for x in all_strengths(filter_part_names(p_data, ['P1'])) if x != None]
	ax.hist(cur_data, bins=prom_bin_list, color=[0.9,0.6,0.6], edgecolor=[0.9,0.6,0.6])
	plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)
	if np.median(cur_data) > 0.0:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	else:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	plt.text(0.02, 0.95, plot_stats, fontsize=14,
 		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
	plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.set_xlim([x_min,x_max])
	ax.set_ylim([0,70])
	# P2
	ax = plt.subplot(3,1,2)
	simpleaxis(ax)
	cur_data = [x for x in all_strengths(filter_part_names(p_data, ['P2'])) if x != None]
	ax.hist(cur_data, bins=prom_bin_list, color=[0.6,0.9,0.6], edgecolor=[0.6,0.9,0.6])
	plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)
	if np.median(cur_data) > 0.0:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	else:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	plt.text(0.02, 0.95, plot_stats, fontsize=14,
 		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
	plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.set_xlim([x_min,x_max])
	ax.set_ylim([0,25])
	# P3
	ax = plt.subplot(3,1,3)
	simpleaxis(ax)
	cur_data = [x for x in all_strengths(filter_part_names(p_data, ['P3'])) if x != None]
	ax.hist(cur_data, bins=prom_bin_list, color=[0.6,0.6,0.9], edgecolor=[0.6,0.6,0.9])
	plt.axvline(x=np.median(cur_data), color='r', linewidth=2.0)
	if np.median(cur_data) > 0.0:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	else:
		plot_stats = 'N = '+str(len(cur_data)) + '\nM = '+str(np.median(cur_data))
	plt.text(0.02, 0.95, plot_stats, fontsize=14,
 		horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
	plt.axvline(x=0, color=(0.2,0.2,0.2), linestyle='--', linewidth=1.5)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.set_xlim([x_min,x_max])
	ax.set_ylim([0,25])
	plt.tight_layout()
	plt.savefig(RESULTS_DIR+'prom_All.pdf')

###############################################################################
# 2. IN/OUT FLUXES
###############################################################################

p_data_parts = split_on_part(p_data)

fig = plt.figure(figsize=(3,6))
ax = plt.subplot(3,1,1)
simpleaxis(ax)
cur_flux_in = [x for x in all_in_fluxes(p_data_parts['P1']) if x != None]
cur_flux_out = np.array([x for x in all_out_fluxes(p_data_parts['P1']) if x != None]) - np.array(cur_flux_in)
ax.scatter(cur_flux_in,cur_flux_out, color=[0.9,0.6,0.6])

ax = plt.subplot(3,1,2)
simpleaxis(ax)
cur_flux_in = [x for x in all_in_fluxes(p_data_parts['P2']) if x != None]
cur_flux_out = np.array([x for x in all_out_fluxes(p_data_parts['P2']) if x != None]) - np.array(cur_flux_in)
ax.scatter(cur_flux_in,cur_flux_out, color=[0.6,0.9,0.6])

ax = plt.subplot(3,1,3)
simpleaxis(ax)
cur_flux_in = [x for x in all_in_fluxes(p_data_parts['P3']) if x != None]
cur_flux_out = np.array([x for x in all_out_fluxes(p_data_parts['P3']) if x != None]) - np.array(cur_flux_in)
ax.scatter(cur_flux_in,cur_flux_out, color=[0.6,0.6,0.9])

plt.tight_layout()
plt.savefig(RESULTS_DIR+'flux_in_out_Promoters.pdf')

# Terminator fluxes
t_data_parts = split_on_part(t_data)

fig = plt.figure(figsize=(3,3))
ax = plt.subplot(1,1,1)
simpleaxis(ax)
cur_flux_in = [x for x in all_in_fluxes(t_data) if x != None]
cur_flux_out = [x for x in all_out_fluxes(t_data) if x != None]
ax.scatter(cur_flux_in,cur_flux_out, s=4.0, color=[0.5,0.5,0.5])
ax.set_xlim([1000,300000])
ax.set_ylim([0.1,300000])
ax.plot(np.arange(0.1, 300000), np.arange(0.1, 300000), color=(0,0,0), linestyle='--', linewidth=1.2)
ax.set_yscale('log')
ax.set_xscale('log')

plt.tight_layout()
plt.savefig(RESULTS_DIR+'flux_in_out_Terminators.pdf')

###############################################################################
# 3. SEQUENCE ANALYSIS
###############################################################################

SEQ_OFFSET = 50

# Extract 20bp up/down stream of every part (act as background set)
p1_data = filter_part_names(p_data, ['P1'])
p2_data = filter_part_names(p_data, ['P2'])
p3_data = filter_part_names(p_data, ['P3'])

p1_seq_data = []
for p in p1_data:
	cur_variant = p[0]
	cur_part_idx = p[1]
	cur_strength = p[-1]
	cur_seq = nifs.extract_seq_range(cur_variant, cur_part_idx, SEQ_OFFSET, SEQ_OFFSET)
	cur_us_seq = cur_seq[0:SEQ_OFFSET]
	cur_ds_seq = cur_seq[-SEQ_OFFSET:]
	p1_seq_data.append([cur_us_seq, cur_ds_seq, cur_seq] + cur_strength)
p1_seq_data = sorted(p1_seq_data, key=lambda x: x[5])

p2_seq_data = []
for p in p2_data:
	cur_variant = p[0]
	cur_part_idx = p[1]
	cur_strength = p[-1]
	cur_seq = nifs.extract_seq_range(cur_variant, cur_part_idx, SEQ_OFFSET, SEQ_OFFSET)
	cur_us_seq = cur_seq[0:SEQ_OFFSET]
	cur_ds_seq = cur_seq[-SEQ_OFFSET:]
	p2_seq_data.append([cur_us_seq, cur_ds_seq, cur_seq] + cur_strength)
p2_seq_data = sorted(p2_seq_data, key=lambda x: x[5])

p3_seq_data = []
for p in p3_data:
	cur_variant = p[0]
	cur_part_idx = p[1]
	cur_strength = p[-1]
	cur_seq = nifs.extract_seq_range(cur_variant, cur_part_idx, SEQ_OFFSET, SEQ_OFFSET)
	cur_us_seq = cur_seq[0:SEQ_OFFSET]
	cur_ds_seq = cur_seq[-SEQ_OFFSET:]
	p3_seq_data.append([cur_us_seq, cur_ds_seq, cur_seq] + cur_strength)
p3_seq_data = sorted(p3_seq_data, key=lambda x: x[5])

t_seq_data = []
for t in t_data:
	cur_variant = t[0]
	cur_part_idx = t[1]
	cur_strength = t[-1]
	cur_seq = nifs.extract_seq_range(cur_variant, cur_part_idx, SEQ_OFFSET, SEQ_OFFSET)
	cur_us_seq = cur_seq[0:SEQ_OFFSET]
	cur_ds_seq = cur_seq[-SEQ_OFFSET:]
	t_seq_data.append([cur_us_seq, cur_ds_seq, cur_seq] + cur_strength)
t_seq_data = sorted(t_seq_data, key=lambda x: x[5])

# Look at highest 10% and lowest 10% of strengths for each Promoter/Terminator
p1_len_10_per = int(len(p1_seq_data) * 0.1)
p1_top_10_per = p1_seq_data[-p1_len_10_per:]
p1_bot_10_per = p1_seq_data[0:p1_len_10_per]

p2_len_10_per = int(len(p2_seq_data) * 0.1)
p2_top_10_per = p2_seq_data[-p2_len_10_per:]
p2_bot_10_per = p2_seq_data[0:p2_len_10_per]

p3_len_10_per = int(len(p3_seq_data) * 0.1)
p3_top_10_per = p3_seq_data[-p3_len_10_per:]
p3_bot_10_per = p3_seq_data[0:p3_len_10_per]

t_len_10_per = int(len(t_seq_data) * 0.1)
t_top_10_per = t_seq_data[-t_len_10_per:]
t_bot_10_per = t_seq_data[0:t_len_10_per]

# Save the reference and top/bottom sets to FASTA file
def save_seqs_to_file (seq_list, seq_idx, name_prefix, out_filename):
	f = open(out_filename, 'w')
	i = 0
	for s in seq_list:
		f.write(str(s[seq_idx]).upper() + '\n')
		i += 1
	f.close()

# Reference sets (all)
save_seqs_to_file(p1_seq_data, 0, 'ref_P1_us_', SEQ_RESULTS_DIR+'ref_p1_us.txt')
save_seqs_to_file(p1_seq_data, 1, 'ref_P1_ds_', SEQ_RESULTS_DIR+'ref_p1_ds.txt')
save_seqs_to_file(p2_seq_data, 0, 'ref_P2_us_', SEQ_RESULTS_DIR+'ref_p2_us.txt')
save_seqs_to_file(p2_seq_data, 1, 'ref_P2_ds_', SEQ_RESULTS_DIR+'ref_p2_ds.txt')
save_seqs_to_file(p3_seq_data, 0, 'ref_P3_us_', SEQ_RESULTS_DIR+'ref_p3_us.txt')
save_seqs_to_file(p3_seq_data, 1, 'ref_P3_ds_', SEQ_RESULTS_DIR+'ref_p3_ds.txt')
save_seqs_to_file(t_seq_data,  0, 'ref_T_us_', SEQ_RESULTS_DIR+'ref_t_us.txt')
save_seqs_to_file(t_seq_data,  1, 'ref_T_ds_', SEQ_RESULTS_DIR+'ref_t_ds.txt')

# Top/Bottom sets
save_seqs_to_file(p1_top_10_per, 0, 'top_P1_us_', SEQ_RESULTS_DIR+'top_p1_us.txt')
save_seqs_to_file(p1_top_10_per, 1, 'top_P1_ds_', SEQ_RESULTS_DIR+'top_p1_ds.txt')
save_seqs_to_file(p1_bot_10_per, 0, 'bot_P1_us_', SEQ_RESULTS_DIR+'bot_p1_us.txt')
save_seqs_to_file(p1_bot_10_per, 1, 'bot_P1_ds_', SEQ_RESULTS_DIR+'bot_p1_ds.txt')

save_seqs_to_file(p2_top_10_per, 0, 'top_P2_us_', SEQ_RESULTS_DIR+'top_p2_us.txt')
save_seqs_to_file(p2_top_10_per, 1, 'top_P2_ds_', SEQ_RESULTS_DIR+'top_p2_ds.txt')
save_seqs_to_file(p2_bot_10_per, 0, 'bot_P2_us_', SEQ_RESULTS_DIR+'bot_p2_us.txt')
save_seqs_to_file(p2_bot_10_per, 1, 'bot_P2_ds_', SEQ_RESULTS_DIR+'bot_p2_ds.txt')

save_seqs_to_file(p3_top_10_per, 0, 'top_P3_us_', SEQ_RESULTS_DIR+'top_p3_us.txt')
save_seqs_to_file(p3_top_10_per, 1, 'top_P3_ds_', SEQ_RESULTS_DIR+'top_p3_ds.txt')
save_seqs_to_file(p3_bot_10_per, 0, 'bot_P3_us_', SEQ_RESULTS_DIR+'bot_p3_us.txt')
save_seqs_to_file(p3_bot_10_per, 1, 'bot_P3_ds_', SEQ_RESULTS_DIR+'bot_p3_ds.txt')

save_seqs_to_file(t_top_10_per, 0, 'top_T_us_', SEQ_RESULTS_DIR+'top_t_us.txt')
save_seqs_to_file(t_top_10_per, 1, 'top_T_ds_', SEQ_RESULTS_DIR+'top_t_ds.txt')
save_seqs_to_file(t_bot_10_per, 0, 'bot_T_us_', SEQ_RESULTS_DIR+'bot_t_us.txt')
save_seqs_to_file(t_bot_10_per, 1, 'bot_T_ds_', SEQ_RESULTS_DIR+'bot_t_ds.txt')


###############################################################################
# CLEAN-UP
###############################################################################

# Clear the plotting cache
plt.close('all')








"""
fig = plt.figure(figsize=(5,5))
ax1 = plt.subplot(1,1,1)
ax1.scatter(all_in_fluxes(p_data), all_strengths(p_data))
ax1.set_title('Promoter [In Flux] vs [Strength]', fontsize=12)



fig = plt.figure(figsize=(5,5))
ax1 = plt.subplot(1,1,1)
ax1.scatter(all_in_fluxes(t_data), all_strengths(t_data))
ax1.set_title('Terminator [In Flux] vs [Strength]', fontsize=12)
"""
"""
print 'TERMINATORS...'
t_data_contexts = split_on_context(t_data)
for c in t_data_contexts.keys():
	cur_data = [x for x in all_strengths(t_data_contexts[c]) if x != None]
	print c, len(cur_data), np.median(cur_data), np.std(cur_data)

print 'PROMOTERS...'
p_data_contexts = split_on_context(filter_part_names(p_data, ['P3', 'P1', 'P2']))
for c in p_data_contexts.keys():
	cur_data = [x for x in all_strengths(p_data_contexts[c]) if x != None]
	print c, len(cur_data), np.median(cur_data), np.std(cur_data)
"""

"""
fig = plt.figure(figsize=(5,5))
ax1 = plt.subplot(1,1,1)
ax1.hist([x for x in all_strengths(p_data) if x != None])
ax1.set_title('p_data [Strength]', fontsize=12)
"""
