#!/usr/bin/env python
"""
Repressor analysis (E.coli DH10B)
"""

import motif_find as mf

import numpy as np
import math
import scipy.stats as stats
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT\n\
			   Jing Zhang <jgzhang@mit.edu>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

###############################################################################
# Parameters
genome_filename = '../data/genomes/DH10B.fasta'
tss_filename = '../data/tss/DH10B_TSS.tsv'
fpkm_filename = '../data/exp/fpkm_data.txt'
metadata_filename = '../data/genomes/Ecoli_metadata.tsv'
results_prefix = '../results/'
repressors = ['AmeR', 'McbR', 'PsrA', 'SrpR', 
              'AmtR', 'ButR', 'QacR', 'TarA',
              'ScbR', 'Orf2',
              'LitR']
#repressors = ['AmtR', 'LitR'] #, 'SrpR']
us_window = 100
ds_window = 100

thresholds = [0.0, 2.5, 5.0]

fmt_axis_outline_width = 1.5
fmt_axis_label_fontsize = 10
###############################################################################

# Load data files
genome_seq = mf.load_genome_seq(genome_filename)
background = mf.genome_background(genome_seq)
genome_len = len(str(genome_seq['DH10B']))
tss = mf.load_tss(tss_filename)
metadata = mf.load_metadata(metadata_filename)
hits_all = {}
hits_tss = {}
motifs = {}
hits_at_tss = {}

# Process each repressor and load all data
for r in repressors:
	print("Processing "+r)
	# Load data
	motifs[r] = mf.load_motifs('../data/motifs/'+r+'.txt')
	hits_all[r] = mf.load_score_hits (results_prefix+r+'_score_hits.csv')
	# Because a repressor has a width that means even further upstream binding sites will play a role
	us_window_corr = us_window+len(motifs[r])
	hits_tss[r] = mf.extract_tss_hits(motifs[r], genome_len, tss, hits_all[r], us_window_corr, ds_window)
	hits_at_tss[r] = mf.hit_per_tss(tss, hits_all[r], motifs[r], genome_len, us_window, ds_window)

def format_axes (ax):
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(fmt_axis_label_fontsize)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(fmt_axis_label_fontsize)

###############################################################################
# SAVE GENERAL STATS

f_out = open(results_prefix+'/_plots/stats.csv', 'w')
f_out.write('repressor,#all_hits,#tss_hits,#tss_hits/#all_hits,#all_hits,#tss_hits,#tss_hits/#all_hits,#all_hits,#tss_hits,#tss_hits/#all_hits\n')
for r in repressors:

	all_hits_weak = mf.filter_hits(hits_all[r], 0.0)
	all_hits_mid = mf.filter_hits(hits_all[r], 2.5)
	all_hits_strong = mf.filter_hits(hits_all[r], 5.0)

	tss_hits_weak = mf.filter_hits(hits_tss[r], 0.0)
	tss_hits_mid = mf.filter_hits(hits_tss[r], 2.5)
	tss_hits_strong = mf.filter_hits(hits_tss[r], 5.0)

	out_data = []
	out_data.append(r)

	out_data.append(str(len(all_hits_weak['DH10B'][:,1])))
	out_data.append(str(len(tss_hits_weak['DH10B'][:,1])))
	out_data.append(str(float(len(tss_hits_weak['DH10B'][:,1]))/len(all_hits_weak['DH10B'][:,1])))

	out_data.append(str(len(all_hits_mid['DH10B'][:,1])))
	out_data.append(str(len(tss_hits_mid['DH10B'][:,1])))
	out_data.append(str(float(len(tss_hits_mid['DH10B'][:,1]))/len(all_hits_mid['DH10B'][:,1])))

	out_data.append(str(len(all_hits_strong['DH10B'][:,1])))
	out_data.append(str(len(tss_hits_strong['DH10B'][:,1])))
	out_data.append(str(float(len(tss_hits_strong['DH10B'][:,1]))/len(all_hits_strong['DH10B'][:,1])))

	f_out.write(','.join(out_data)+'\n')
f_out.close()

###############################################################################
# Essential genes that are hit

ess_locus_ids = []
for gene in metadata.keys():
	if metadata[gene][-1] == 'Y':
		for locus in tss['DH10B'].keys():
			if gene in tss['DH10B'][locus][0]:
				ess_locus_ids.append(locus)

ess_bindings = {}
for r in repressors:
	# check to see how many of the strong binding sites are essential
	ess_bindings[r] = []
	for ess_locus in ess_locus_ids:
		if ess_locus in hits_at_tss[r]['DH10B'].keys():
			ess_bindings[r].append([ess_locus] + hits_at_tss[r]['DH10B'][ess_locus])








###############################################################################
# Changes in FPKM values

a = False
if a == True:

	fpkms = mf.load_fpkms(fpkm_filename)
	# Go through each repressor extract expression and expression where no hit and hit
	rep_corr_hits = {}
	rep_corr_nohits = {}
	for r in repressors:
		rep_corr_hits[r] = []
		rep_corr_nohits[r] = []
		rep_locus_id = 'SYNTHETIC_'+r
		x = fpkms['DH10B'][rep_locus_id]
		
		for c in fpkms.keys():
			for locus in tss[c].keys():
				if locus != rep_locus_id and locus[0] != 'S':
					hit_num = 0
					if locus in hits_at_tss[r][c].keys():
						for el in hits_at_tss[r][c][locus]:
							if el[2] > 5.0:
								hit_num += 1

					if locus in fpkms[c].keys():
						y = fpkms[c][locus]
						rho, p_val = stats.spearmanr(x,y)
						if not math.isnan(rho):
							if hit_num > 0:
								rep_corr_hits[r].append([locus, rho, p_val])
							else:
								rep_corr_nohits[r].append([locus, rho, p_val])

	#print 'Repressor: '+r+', no hit rho = '+str(np.mean(rep_corr_nohits[r]))+', hit rho = '+str(np.mean(rep_corr_hits[r]))
	#print 'Repressor: '+r+', no hit rho = '+str(np.median(rep_corr_nohits[r]))+', hit rho = '+str(np.median(rep_corr_hits[r]))

###############################################################################
# DISTRIBUTIONS OF HIT NUMBERS AND SCORES

col_dist_all = (0.2,0.2,0.8,0.4)
col_dist_tss = (0.8,0.2,0.2,0.4)

gs = gridspec.GridSpec(11, 2, width_ratios=[1,2])
fig = plt.figure(figsize=(6.5,12))
rep_idx = 0
for r in repressors:

	ax = plt.subplot(gs[(rep_idx*2)])
	cur_threshold = 0.0
	x_vals, y_vals = mf.binding_trace(tss, hits_tss[r], motifs[r], genome_len, us_window, ds_window, cur_threshold)
	ax.plot(x_vals, y_vals, color=(0,0,0), linewidth=1.0, zorder=70)
	y_max = max(y_vals)
	#cur_threshold = 1.5
	#x_vals, y_vals = mf.binding_trace(tss, hits_tss[r], motifs[r], genome_len, us_window, ds_window, cur_threshold)
	#ax.plot(x_vals, y_vals, color=(0.2,0.2,0.2), linewidth=1.0, zorder=80)
	cur_threshold = 3.0
	x_vals, y_vals = mf.binding_trace(tss, hits_tss[r], motifs[r], genome_len, us_window, ds_window, cur_threshold)
	ax.plot(x_vals, y_vals, color=(0.4,0.4,0.4), linewidth=1.0, zorder=90)
	
	# Format data
	plt.axvline(x=0, linestyle='-', linewidth=1.5, color=(0,0,0))
	#plt.axhline(y=0, linestyle='-', linewidth=1.0, color=(0,0,0))
	ax.fill_between([-55,20], -999, 999, linewidth=0, facecolor=(0.85,0.85,0.85), zorder=0)
	ax.set_ylim([-(y_max*0.1), y_max*1.1])
	ax.set_xlim([-us_window, ds_window])
	ax.set_xticks([-us_window, 0, ds_window])
	if rep_idx < len(repressors)-1:
		ax.xaxis.set_ticklabels([])
	plt.tight_layout()
	format_axes(ax)
	plt.subplots_adjust(hspace=.2)

	###############################################################################

	ax = plt.subplot(gs[(rep_idx*2)+1])
	# Make sure histograms are comparable (density functions)
	bin_list = np.linspace(0, 20, 50)
	ax.hist(hits_all[r]['DH10B'][:,1], bins=bin_list, normed=True, linewidth=1.25, color=(0,0,0), facecolor=(0.0,0.0,0.0,0.0), histtype='stepfilled', zorder=-1)
	if len(hits_tss[r]['DH10B']) != 0:
		ax.hist(hits_tss[r]['DH10B'][:,1], bins=bin_list, normed=True, linewidth=0, facecolor=(0.6,0.6,0.6), alpha=1.0, histtype='stepfilled', zorder=-2)
	if rep_idx < len(repressors)-1:
		ax.xaxis.set_ticklabels([])
	ax.set_xlim([0,20])
	
	if r == 'AmtR':
		ax.set_ylim([0,2.5])
	else:
		ax.set_ylim([0,0.7])
	plt.tight_layout()
	format_axes(ax)
	plt.subplots_adjust(hspace=.2)
	rep_idx += 1

fig.savefig(results_prefix+'/_plots/hit_dists.pdf', transparent=True)
plt.close('all')


