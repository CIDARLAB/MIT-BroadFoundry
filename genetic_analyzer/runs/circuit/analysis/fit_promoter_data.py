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
OUTPUT_PREFIX = './plots_fit/'

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

# Colour maps to use for the genes
cmap = {}
cmap['AmtR'] = (1.0,0.75,0.17) # 255, 193, 43
cmap['LitR'] = (0.38,0.82,0.32) # 98, 209, 83
cmap['BM3R1'] = (0.95,0.30,0.25) # 242, 78, 65
cmap['SrpR'] = (0.38,0.65,0.87) # 97, 165, 223
cmap['PhlF'] = (0.55,0.35,0.64) # 141, 89, 163
cmap['YFP'] = (0.98,0.97,0.35) # 250, 248, 89

reads_to_YFP_fac = 7.18228937996
fac_reduce_size = 1000.0/reads_to_YFP_fac

###############################################################################
# PERFORM FITTING
###############################################################################

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
		# Numbers are generally in read-depths (reduce to reduce squared values being too large)
		err_diffs.append((exp_val-fit_val)/fac_reduce_size)
	# Return SSE
	return np.sum(np.power(err_diffs, 2.0))


def induced_err_func (x, exp_data, samples_to_include, chrom_to_fit, PU_to_fit, P_input_states, P_names):		
	# Extract the parameters into more user friendly dict
	fit_params = {}
	for p_idx in range(len(P_names)):
		p_name = P_names[p_idx]
		fit_params[p_name] = {}
		fit_params[p_name]['Pmin'] = x[p_idx*2]
		fit_params[p_name]['Pmin_inc'] = x[(p_idx*2)+1]
	# Only calculate error based on samples specified
	err_diffs = []
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		exp_val = exp_data[sample][chrom_to_fit][PU_to_fit]
		# Calculate fitted output assuming additive fluxes
		fit_outs = []
		off_fac = 1.0
		for p_idx in range(len(P_names)):
			is_on = False
			if P_input_states[p_idx][sample] == 1:
				is_on = True
			cur_p = P_names[p_idx]
			if is_on == True:
				fit_outs.append(fit_params[cur_p]['Pmin']+fit_params[cur_p]['Pmin_inc'])
			else:
				off_fac = 1000.0
				fit_outs.append(fit_params[cur_p]['Pmin'])
		fit_val = np.sum(fit_outs)
		# Numbers are generally in read-depths (reduce to reduce squared values being too large)
		err_diffs.append(((exp_val-fit_val)*off_fac)/fac_reduce_size)
	# Return SSE
	return np.sum(np.power(err_diffs, 2.0))

def fit_induced_to_data (exp_data, samples_to_include, chrom_to_fit, PU_to_fit, P_input_states, P_names):
	num_of_promoters = len(P_names)
	# Parameters for fit (Pmin, Pmax)
	x0 = np.zeros(num_of_promoters*2)
	# Some contraints (keep everything positive)
	bnds = []
	for p in range(num_of_promoters):
		# Initial conditions (start realistic to improve fitting)
		x0[p*2] = 0.0
		x0[(p*2)+1] = 1000.0
		# Set bounds
		bnds.append((0.0, None))
		bnds.append((0.0, None))
	# methods = BFGS, nelder-mead, Powell, TNC, SLSQP,  L-BFGS-B
	res = scipy.optimize.minimize(induced_err_func, x0, args=(exp_data, samples_to_include, 
		                          chrom_to_fit, PU_to_fit, P_input_states, P_names),
								  bounds=bnds,
		                          method='SLSQP', jac=False,
		                          options={'disp': True, 'maxiter': 10000})
	# Reformat results into dictionary
	params = {}
	for p_idx in range(len(P_names)):
		p_name = P_names[p_idx]
		params[p_name] = {}
		params[p_name]['Pmin'] = res.x[p_idx*2]
		params[p_name]['Pmin_inc'] = res.x[(p_idx*2)+1]
	return params, res.success


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
		bnds.append((1.0, 3.9))
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
			cur_strength = float(row[5]) * reads_to_YFP_fac
			if cur_sample not in exp_data.keys():
				exp_data[cur_sample] = {}
			if cur_chrom not in exp_data[cur_sample].keys():
				exp_data[cur_sample][cur_chrom] = {}
			if cur_strength < 0.0:
				exp_data[cur_sample][cur_chrom][cur_part] = 0.0
			else:
				exp_data[cur_sample][cur_chrom][cur_part] = cur_strength
	return exp_data


def analyse_promoters (exp_data, samples_to_include):
	ind_states = {}
	ind_states['pTac'] = [0,1,0,1,0,1,0,1]
	ind_states['pTet'] = [0,0,1,1,0,0,1,1]
	ind_states['pBAD'] = [0,0,0,0,1,1,1,1]

	pTac_ind_states = {'tube_1': 0, 'tube_2': 1, 'tube_3': 0, 'tube_4': 1, 'tube_5': 0, 'tube_6': 1, 'tube_7': 0, 'tube_8': 1,
	                   'flask_1': 0, 'flask_2': 1, 'flask_3': 0, 'flask_4': 1, 'flask_5': 0, 'flask_6': 1, 'flask_7': 0, 'flask_8': 1}

	pTet_ind_states = {'tube_1': 0, 'tube_2': 0, 'tube_3': 1, 'tube_4': 1, 'tube_5': 0, 'tube_6': 0, 'tube_7': 1, 'tube_8': 1, 
	                   'flask_1': 0, 'flask_2': 0, 'flask_3': 1, 'flask_4': 1, 'flask_5': 0, 'flask_6': 0, 'flask_7': 1, 'flask_8': 1}

	pBAD_ind_states = {'tube_1': 0, 'tube_2': 0, 'tube_3': 0, 'tube_4': 0, 'tube_5': 1, 'tube_6': 1, 'tube_7': 1, 'tube_8': 1,
	                   'flask_1': 0, 'flask_2': 0, 'flask_3': 0, 'flask_4': 0, 'flask_5': 1, 'flask_6': 1, 'flask_7': 1, 'flask_8': 1}

	pu_Tac_Tet1_params, pu_Tac_Tet1_suc = fit_induced_to_data (exp_data, samples_to_include,
		'0x58v50', 'PU-pTac-pTet1', 
		[pTac_ind_states, pTet_ind_states],
		['pTac', 'pTet1'])

	pu_BAD1_Tet2_params, pu_Tac_Tet1_suc = fit_induced_to_data (exp_data, samples_to_include,
		'0x58v50', 'PU-pBAD1-pTet2', 
		[pBAD_ind_states, pTet_ind_states],
		['pBAD1', 'pTet2'])

	pu_BAD2_params, pu_Tac_Tet1_suc = fit_induced_to_data (exp_data, samples_to_include,
		'0x58v50', 'PU-pBAD2', 
		[pBAD_ind_states],
		['pBAD2'])

	pu_BM3R1_AmtR_params, pu_BM3R1_AmtR_suc = fit_repressors_to_data(exp_data, samples_to_include, 
		'0x58v50', 'PU-pBM3R1-pAmtR', 
		['0x58v50', '0x58v50'], ['PU-pBAD2', 'PU-pTac-pTet1'], 
		['pBM3R1', 'pAmtR'])

	pu_SrpR_LitR_params, pu_SrpR_LitR_suc = fit_repressors_to_data(exp_data, samples_to_include, 
		'0x58v50', 'PU-pSrpR-pLitR', 
		['0x58v50', '0x58v50'], ['PU-pBM3R1-pAmtR', 'PU-pBAD1-pTet2'], 
		['pSrpR', 'pLitR'])

	pu_PhlF_params, pu_BM3R1_AmtR_suc = fit_repressors_to_data(exp_data, samples_to_include, 
		'0x58v50', 'PU-pPhlF', 
		['0x58v50'], ['PU-pSrpR-pLitR'], 
		['pPhlF'])

	combined_params = {}
	for k in pu_BM3R1_AmtR_params.keys():
		combined_params[k] = pu_BM3R1_AmtR_params[k]
	for k in pu_SrpR_LitR_params.keys():
		combined_params[k] = pu_SrpR_LitR_params[k]
	for k in pu_PhlF_params.keys():
		combined_params[k] = pu_PhlF_params[k]

	for k in pu_Tac_Tet1_params.keys():
		combined_params[k] = pu_Tac_Tet1_params[k]
	for k in pu_BAD1_Tet2_params.keys():
		combined_params[k] = pu_BAD1_Tet2_params[k]
	for k in pu_BAD2_params.keys():
		combined_params[k] = pu_BAD2_params[k]

	ind_data = {}
	pred_out1 = 0.0
	pred_out2 = 0.0
	ind_data['Tac'] = [[],[]]
	ind_data['Tet1'] = [[],[]]
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		# Calculate inputs
		in1 = pTac_ind_states[sample]
		in2 = pTet_ind_states[sample]
		ind_data['Tac'][0].append(in1)
		ind_data['Tet1'][0].append(in2)
		# Calculate expected outputs
		if pTac_ind_states[sample] == 1:
			pred_out1 = combined_params['pTac']['Pmin']+combined_params['pTac']['Pmin_inc']
		else:
			pred_out1 = combined_params['pTac']['Pmin']
		if pTet_ind_states[sample] == 1:
			pred_out2 = combined_params['pTet1']['Pmin']+combined_params['pTet1']['Pmin_inc']
		else:
			pred_out2 = combined_params['pTet1']['Pmin']
		# Fractional contribution of each promoter
		fac1 = pred_out1/(pred_out1+pred_out2)
		fac2 = pred_out2/(pred_out1+pred_out2)
		# Split output by predicted contribution factors
		out1 = fac1 * exp_data[sample]['0x58v50']['PU-pTac-pTet1']
		out2 = fac2 * exp_data[sample]['0x58v50']['PU-pTac-pTet1']
		ind_data['Tac'][1].append(out1)
		ind_data['Tet1'][1].append(out2)


	ind_data['BAD1'] = [[],[]]
	ind_data['Tet2'] = [[],[]]
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		# Calculate inputs
		in1 = pBAD_ind_states[sample]
		in2 = pTet_ind_states[sample]
		ind_data['BAD1'][0].append(in1)
		ind_data['Tet2'][0].append(in2)
		# Calculate expected outputs
		if pBAD_ind_states[sample] == 1:
			pred_out1 = combined_params['pBAD1']['Pmin']+combined_params['pBAD1']['Pmin_inc']
		else:
			pred_out1 = combined_params['pBAD1']['Pmin']
		if pTet_ind_states[sample] == 1:
			pred_out2 = combined_params['pTet2']['Pmin']+combined_params['pTet2']['Pmin_inc']
		else:
			pred_out2 = combined_params['pTet2']['Pmin']
		# Fractional contribution of each promoter
		fac1 = pred_out1/(pred_out1+pred_out2)
		fac2 = pred_out2/(pred_out1+pred_out2)
		# Split output by predicted contribution factors
		out1 = fac1 * exp_data[sample]['0x58v50']['PU-pBAD1-pTet2']
		out2 = fac2 * exp_data[sample]['0x58v50']['PU-pBAD1-pTet2']
		ind_data['BAD1'][1].append(out1)
		ind_data['Tet2'][1].append(out2)

	ind_data['BAD2'] = [[],[]]
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		# Calculate inputs
		in1 = pBAD_ind_states[sample]
		ind_data['BAD2'][0].append(in1)
		# Calculate expected outputs
		if pBAD_ind_states[sample] == 1:
			pred_out1 = combined_params['pBAD2']['Pmin']+combined_params['pBAD2']['Pmin_inc']
		else:
			pred_out1 = combined_params['pBAD2']['Pmin']
		# Split output by predicted contribution factors
		out1 = exp_data[sample]['0x58v50']['PU-pBAD2']
		ind_data['BAD2'][1].append(out1)

	rep_data = {}
	rep_data['BM3R1'] = [[],[]]
	rep_data['AmtR'] = [[],[]]
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		# Calculate inputs
		in1 = exp_data[sample]['0x58v50']['PU-pBAD2']
		in2 = exp_data[sample]['0x58v50']['PU-pTac-pTet1']
		rep_data['BM3R1'][0].append(in1)
		rep_data['AmtR'][0].append(in2)
		# Calculate expected outputs
		pred_out1 = hill_func(in1, 
			                  pu_BM3R1_AmtR_params['pBM3R1']['Pmin'], 
			                  pu_BM3R1_AmtR_params['pBM3R1']['Pmin_inc'], 
			                  pu_BM3R1_AmtR_params['pBM3R1']['K'],
			                  pu_BM3R1_AmtR_params['pBM3R1']['n'],
			                  repress=True)
		pred_out2 = hill_func(in2, 
			                  pu_BM3R1_AmtR_params['pAmtR']['Pmin'], 
			                  pu_BM3R1_AmtR_params['pAmtR']['Pmin_inc'], 
			                  pu_BM3R1_AmtR_params['pAmtR']['K'],
			                  pu_BM3R1_AmtR_params['pAmtR']['n'],
			                  repress=True)
		# Fractional contribution of each promoter
		fac1 = pred_out1/(pred_out1+pred_out2)
		fac2 = pred_out2/(pred_out1+pred_out2)
		# Split output by predicted contribution factors
		out1 = fac1 * exp_data[sample]['0x58v50']['PU-pBM3R1-pAmtR']
		out2 = fac2 * exp_data[sample]['0x58v50']['PU-pBM3R1-pAmtR']
		rep_data['BM3R1'][1].append(out1)
		rep_data['AmtR'][1].append(out2)

	rep_data['SrpR'] = [[],[]]
	rep_data['LitR'] = [[],[]]
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		# Calculate inputs
		in1 = exp_data[sample]['0x58v50']['PU-pBM3R1-pAmtR']
		in2 = exp_data[sample]['0x58v50']['PU-pBAD1-pTet2']
		rep_data['SrpR'][0].append(in1)
		rep_data['LitR'][0].append(in2)
		# Calculate expected outputs
		pred_out1 = hill_func(in1, 
			                  pu_SrpR_LitR_params['pSrpR']['Pmin'], 
			                  pu_SrpR_LitR_params['pSrpR']['Pmin_inc'], 
			                  pu_SrpR_LitR_params['pSrpR']['K'],
			                  pu_SrpR_LitR_params['pSrpR']['n'],
			                  repress=True)
		pred_out2 = hill_func(in2, 
			                  pu_SrpR_LitR_params['pLitR']['Pmin'], 
			                  pu_SrpR_LitR_params['pLitR']['Pmin_inc'], 
			                  pu_SrpR_LitR_params['pLitR']['K'],
			                  pu_SrpR_LitR_params['pLitR']['n'],
			                  repress=True)
		# Fractional contribution of each promoter
		fac1 = pred_out1/(pred_out1+pred_out2)
		fac2 = pred_out2/(pred_out1+pred_out2)
		# Split output by predicted contribution factors
		out1 = fac1 * exp_data[sample]['0x58v50']['PU-pSrpR-pLitR']
		out2 = fac2 * exp_data[sample]['0x58v50']['PU-pSrpR-pLitR']
		rep_data['SrpR'][1].append(out1)
		rep_data['LitR'][1].append(out2)

	rep_data['PhlF'] = [[],[]]
	for sample in samples_to_include:
		sample_data = exp_data[sample]
		# Calculate inputs
		in1 = exp_data[sample]['0x58v50']['PU-pSrpR-pLitR']
		rep_data['PhlF'][0].append(in1)
		# Calculate expected outputs
		pred_out1 = hill_func(in1, 
			                  pu_PhlF_params['pPhlF']['Pmin'], 
			                  pu_PhlF_params['pPhlF']['Pmin_inc'], 
			                  pu_PhlF_params['pPhlF']['K'],
			                  pu_PhlF_params['pPhlF']['n'],
			                  repress=True)
		# Split output by predicted contribution factors
		out1 = exp_data[sample]['0x58v50']['PU-pPhlF']
		rep_data['PhlF'][1].append(out1)

	return combined_params, ind_data, rep_data

###############################################################################
# PLOT RESULTS
###############################################################################

def plot_char_repressor (Pmin, Pmax, K, n, col, x_range, y_range, rep_data, filename_out, marker_shape='o'):
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	# Generate characterisation curve
	#x = np.linspace(x_range[0],x_range[1],1000)
	x = np.logspace(-6,6,500)
	y = []
	for el in x:
		y.append(hill_func(el, Pmin, Pmax, K, n, repress=True))
	y = np.array(y)
	# Add annotations for characterisation
	annot_line_col_1 = (1,1,1)
	annot_line_col_2 = (0,0,0)
	plt.axhline(y=Pmax, color=annot_line_col_2, linewidth=fmt_char_data_linewidth)
	plt.axhline(y=Pmin, color=annot_line_col_2, linewidth=fmt_char_data_linewidth)
	plt.axvline(x=K, color=annot_line_col_2, linestyle='-', linewidth=fmt_char_data_linewidth)
	ax.plot(x,y, linewidth=fmt_char_line_width, color=col)
	
	ax.scatter(rep_data[0], rep_data[1], marker=marker_shape, zorder=100, s=40, color=(1,1,1,0), edgecolor=(0,0,0), linewidth=2.0)
	# Format the axis
	plt.ylim(y_range)
	plt.xlim(x_range)
	ax.set_xscale('symlog', linthreshx=100, linscalex=0.5)
	ax.set_yscale('symlog', linthreshy=100, linscaley=0.5)
	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size, pad=2)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	plt.savefig(filename_out, transparent=True)

def plot_char_inducer (Pmin, Pmax, x_range, y_range, ind_data, filename_out, marker_shape='o'):
	fig = plt.figure(figsize=(2.5,2.5))
	ax = plt.subplot(1,1,1)
	fmt_char_data_linewidth = 1.5
	fmt_axis_outline_width = fmt_edge_width
	fmt_char_line_width = 5.0
	# Add annotations for characterisation
	annot_line_col_1 = (1,1,1)
	annot_line_col_2 = (0,0,0)
	ax.plot([1,4],[Pmin,Pmin], linewidth=fmt_char_line_width, color=(0.5,0.5,0.5))
	ax.plot([6,9],[Pmax,Pmax], linewidth=fmt_char_line_width, color=(0.5,0.5,0.5))
  
	x_vals = np.array(ind_data[0], dtype=np.float)
	x_vals[x_vals==0] = 2.5
	x_vals[x_vals==1] = 7.5
	ax.scatter(x_vals, np.array(ind_data[1]), marker=marker_shape, zorder=100, s=40, color=(1,1,1,0), edgecolor=(0,0,0), linewidth=2.0)
	
	# Format the axis
	ax.set_xlim(x_range)
	ax.set_ylim(y_range)
	ax.set_yscale('symlog', linthreshy=100, linscaley=0.5)

	ax.tick_params(axis='x', labelsize=fmt_label_size)
	ax.tick_params(axis='y', labelsize=fmt_label_size, pad=2)
	ax.get_xaxis().set_ticks([])
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
		ax.spines[axis].set_linewidth(fmt_axis_outline_width)
	plt.subplots_adjust(left=0.18, right=0.90, top=0.90, bottom=0.18)
	plt.savefig(filename_out, transparent=True)

def save_fitted_data (combined_params, filename_out):
	f_out = open(filename_out, 'w')
	for p in combined_params.keys():
		ks = combined_params[p].keys()
		vals = [p]
		for k in ks:
			vals.append(k)
			vals.append(str(combined_params[p][k]))
		f_out.write('\t'.join(vals)+'\n')
	f_out.close()

###############################################################################
# WORKING - TUBE DATA
###############################################################################

exp_data = load_promoter_perf_per_sample('promoter.profile.perf.txt')

samples_to_include = []
for i in range(1,9): #[1,2,3,4,5,6]: #range(1,9):range(1,9):
	samples_to_include.append('tube_'+str(i))

combined_params, ind_data, rep_data = analyse_promoters(exp_data, samples_to_include)

for rep in ['AmtR', 'LitR', 'SrpR', 'BM3R1', 'PhlF']:
	plot_char_repressor (combined_params['p'+rep]['Pmin'],
		                 combined_params['p'+rep]['Pmin']+combined_params['p'+rep]['Pmin_inc'], 
		                 combined_params['p'+rep]['K'],
		                 combined_params['p'+rep]['n'],
		                 cmap[rep], [0.0, 1000000], [-80.0 ,1000000],
		                 rep_data[rep],
		                 OUTPUT_PREFIX+rep+'_tube.pdf', marker_shape='o')

for pro in ['Tac', 'Tet1', 'Tet2', 'BAD1', 'BAD2']:
	plot_char_inducer (combined_params['p'+pro]['Pmin'],
		               combined_params['p'+pro]['Pmin']+combined_params['p'+pro]['Pmin_inc'], 
		               [0.0, 10.0], [-80.0 ,1000000],
		               ind_data[pro],
		               OUTPUT_PREFIX+pro+'_tube.pdf', marker_shape='o')

save_fitted_data (combined_params, OUTPUT_PREFIX+'tube_params.fit.txt')

###############################################################################
# BROKEN - FLASK DATA
###############################################################################

samples_to_include = []
for i in range(1,9):
	samples_to_include.append('flask_'+str(i))

combined_params, ind_data, rep_data = analyse_promoters(exp_data, samples_to_include)

for rep in ['AmtR', 'LitR', 'SrpR', 'BM3R1', 'PhlF']:
	plot_char_repressor (combined_params['p'+rep]['Pmin'],
		                 combined_params['p'+rep]['Pmin']+combined_params['p'+rep]['Pmin_inc'], 
		                 combined_params['p'+rep]['K'],
		                 combined_params['p'+rep]['n'],
		                 cmap[rep], [0.0, 1000000], [-80.0 ,1000000],
		                 rep_data[rep],
		                 OUTPUT_PREFIX+rep+'_flask.pdf', marker_shape='x')

for pro in ['Tac', 'Tet1', 'Tet2', 'BAD1', 'BAD2']:
	plot_char_inducer (combined_params['p'+pro]['Pmin'],
		               combined_params['p'+pro]['Pmin']+combined_params['p'+pro]['Pmin_inc'], 
		               [0.0, 10.0], [-80.0 ,1000000],
		               ind_data[pro],
		               OUTPUT_PREFIX+pro+'_flask.pdf', marker_shape='x')

save_fitted_data (combined_params, OUTPUT_PREFIX+'flask_params.fit.txt')
