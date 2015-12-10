#!/usr/bin/env python
"""
	Fit the data to the promoter units to deconvolve their contributions.
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import csv
import scipy.optimize

matplotlib.rcParams['lines.dash_joinstyle']  = 'miter'
matplotlib.rcParams['lines.dash_capstyle']   = 'butt'
matplotlib.rcParams['lines.solid_joinstyle'] = 'miter'
matplotlib.rcParams['lines.solid_capstyle']  = 'projecting'
# Make text editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype']          = 42

DATA_PREFIX = '../data'
RESULTS_PREFIX = '../results'
OUTPUT_PREFIX = './plots_profile'

MIN_READS_TERM = 500.0
MIN_READS_RIBO = 500.0

tube_samples = ['tube_1', 'tube_2', 'tube_3', 'tube_4', 'tube_5', 'tube_6', 'tube_7', 'tube_8']
flask_samples = ['flask_1', 'flask_2', 'flask_3', 'flask_4', 'flask_5', 'flask_6', 'flask_7', 'flask_8']

plt.rcParams['ytick.major.pad']='1' # 5 for all but repressor graphs
plt.rcParams['xtick.major.pad']='5'
fmt_label_size = 13.5
fmt_edge_width = 2.5
trace_height = 2.5

# Inducer states
ind_states = {}
ind_states['pTac'] = [0,1,0,1,0,1,0,1]
ind_states['pTet'] = [0,0,1,1,0,0,1,1]
ind_states['pBAD'] = [0,0,0,0,1,1,1,1]


def plot_promoter_profile (profiles, trace_region, part_region, filename_out):
	loc = ticker.MultipleLocator(base=300.0)
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	annot_line_col = (0,0,0)
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	for s in profiles:
		profile_data = extract_profile_region(profiles[s], trace_region[0], trace_region[1], trace_region[2])
		plt.axvspan(part_region[1]-trace_region[1], part_region[2]-trace_region[1], facecolor=(0.8,0.8,0.8), linewidth=0, zorder=-10)
		ax.plot(range(trace_region[2]-trace_region[1]), np.array(profile_data[0]), color=(0.0,0.0,0.0), alpha=1, linewidth=1.5)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	ax.set_xlim([0,trace_region[2]-trace_region[1]])
	ax.set_xscale('linear')
	ax.xaxis.set_major_locator(loc)
	ax.set_ylim([0,1000000])
	ax.set_yscale('symlog', linthreshy=10)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')




def hill_func (x, Pmin, Pmin_inc, K, n, repress=False):
	if x < 0.0:
		x = 0.0
	if repress == True:
		return Pmin + (Pmin+Pmin_inc)*( math.pow(K,n) / (math.pow(K,n)+math.pow(x,n)) )
	else: 
		return Pmin + (Pmin+Pmin_inc)*( math.pow(x,n) / (math.pow(K,n)+math.pow(x,n)) )


def repressor_err_func (x, exp_data, samples_to_include, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names):		
	# Extract the parameters into more user friendly dict
	fit_params = {}
	for p_idx in range(len(P_names)):
		p_name = P_names[p_idx]
		fit_params[p_name] = {}
		fit_params[p_name]['Pmin'] = x[p_idx*4]
		fit_params[p_name]['Pmin_inc'] = x[(p_idx*4)+1]
		fit_params[p_name]['K'] = x[(p_idx*4)+2]
		fit_params[p_name]['n'] = x[(p_idx*4)+3]
	# Only calculate error based on samples specified
	err_diffs = []
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		exp_val = exp_data[sample][chrom_to_fit][PU_to_fit]
		# Calculate fitted output assuming additive fluxes
		fit_outs = []
		for p_idx in range(len(P_names)):
			input_val = exp_data[sample][chrom_inputs[p_idx]][PU_inputs[p_idx]]
			cur_p = P_names[p_idx]
			fit_outs.append( hill_func(input_val, 
				                       fit_params[cur_p]['Pmin'], 
				                       fit_params[cur_p]['Pmin_inc'],
				                       fit_params[cur_p]['K'],
				                       fit_params[cur_p]['n'],
				                       repress=True) )
		fit_val = np.sum(fit_outs)
		err_diffs.append((exp_val-fit_val)/1000.0)
	# Return SSE
	return np.sum(np.power(err_diffs, 2.0))


def fit_repressors_to_data (exp_data, samples_to_include, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names):
	num_of_promoters = len(P_names)
	# Parameters for fit (Pmin, Pmax, K, n)
	x0 = np.zeros(num_of_promoters*4)
	# Some contraints (keep everything positive)
	bnds = []
	for p in range(num_of_promoters):
		# Initial conditions (start realistic to improve fitting)
		x0[p*4] = 0.0
		x0[(p*4)+1] = 1000.0
		x0[(p*4)+2] = 100.0
		x0[(p*4)+3] = 2.0
		# Set bounds
		bnds.append((0.0, None))
		bnds.append((0.0, None))
		bnds.append((1.0, None))
		bnds.append((0.1, 5.0))
	# methods = BFGS, nelder-mead, Powell, TNC, SLSQP,  L-BFGS-B
	res = scipy.optimize.minimize(repressor_err_func, x0, args=(exp_data, samples_to_include, 
		                          chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names),
								  bounds=bnds,
		                          method='SLSQP', jac=False,
		                          options={'disp': True, 'maxiter': 10000})
	# Reformat results into dictionary
	params = {}
	for p_idx in range(len(P_names)):
		p_name = P_names[p_idx]
		params[p_name] = {}
		params[p_name]['Pmin'] = res.x[p_idx*4]
		params[p_name]['Pmin_inc'] = res.x[(p_idx*4)+1]
		params[p_name]['K'] = res.x[(p_idx*4)+2]
		params[p_name]['n'] = res.x[(p_idx*4)+3]
	return params, res.success


def load_promoter_perf_per_sample (filename):
	exp_data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	data_reader.next()
	for row in data_reader:
		if len(row) > 5:
			cur_chrom = row[0]
			cur_part = row[1]
			cur_sample = row[2]
			cur_strength = float(row[5])
			if cur_sample not in exp_data.keys():
				exp_data[cur_sample] = {}
			if cur_chrom not in exp_data[cur_sample].keys():
				exp_data[cur_sample][cur_chrom] = {}
			exp_data[cur_sample][cur_chrom][cur_part] = cur_strength
	return exp_data


samples_to_include = []
for i in range(1,9):
	samples_to_include.append('tube_'+str(i))
exp_data = load_promoter_perf_per_sample('promoter.profile.perf.txt')

pu_BM3R1_AmtR_params, pu_BM3R1_AmtR_suc = fit_repressors_to_data(exp_data, samples_to_include, 
	'0x58v50', 'PU-pBM3R1-pAmtR', 
	['0x58v50', '0x58v50'], ['PU-pBAD2', 'PU-pTac-pTet1'], 
	['pBM3R1', 'pAmtR'])
print pu_BM3R1_AmtR_params

pu_SrpR_LitR_params, pu_SrpR_LitR_suc = fit_repressors_to_data(exp_data, samples_to_include, 
	'0x58v50', 'PU-pSrpR-pLitR', 
	['0x58v50', '0x58v50'], ['PU-pBM3R1-pAmtR', 'PU-pBAD1-pTet2'], 
	['pSrpR', 'pLitR'])
print pu_SrpR_LitR_params

pu_PhlF_params, pu_BM3R1_AmtR_suc = fit_repressors_to_data(exp_data, samples_to_include, 
	'0x58v50', 'PU-pPhlF', 
	['0x58v50'], ['PU-pSrpR-pLitR'], 
	['pPhlF'])
print pu_PhlF_params







