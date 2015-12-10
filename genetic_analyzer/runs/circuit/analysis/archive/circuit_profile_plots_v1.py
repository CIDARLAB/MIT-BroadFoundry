#!/usr/bin/env python
"""
	Plot transcription profile plots for paper.
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
import fragmentation_model as fm

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

col_sense_rd = (0.5,0.5,0.5)
col_off_state = (0.95,0.30,0.25)
col_on_state = (0.38,0.65,0.87)

# Colour maps to use for the genes
cmap = {}
cmap['AmtR'] = (1.0,0.75,0.17) # 255, 193, 43
cmap['LitR'] = (0.38,0.82,0.32) # 98, 209, 83
cmap['BM3R1'] = (0.95,0.30,0.25) # 242, 78, 65
cmap['SrpR'] = (0.38,0.65,0.87) # 97, 165, 223
cmap['PhlF'] = (0.55,0.35,0.64) # 141, 89, 163
cmap['YFP'] = (0.98,0.97,0.35) # 250, 248, 89

cmap['pAmtR'] = (1.0,0.75,0.17) # 255, 193, 43
cmap['pLitR'] = (0.38,0.82,0.32) # 98, 209, 83
cmap['pBM3R1'] = (0.95,0.30,0.25) # 242, 78, 65
cmap['pSrpR'] = (0.38,0.65,0.87) # 97, 165, 223
cmap['pPhlF'] = (0.55,0.35,0.64) # 141, 89, 163
cmap['pYFP'] = (0.98,0.97,0.35) # 250, 248, 89

cmap['pTac'] = (0,0,0)
cmap['pTet1'] = (0,0,0)
cmap['pBAD1'] = (0,0,0)
cmap['pTet2'] = (0,0,0)
cmap['pBAD2'] = (0,0,0)

# Terminator colours (all black)
cmap['L3S2P55'] = (0.0,0.0,0.0) # AmtR
cmap['L3S2P24'] = (0.0,0.0,0.0) # LitR
cmap['L3S2P11'] = (0.0,0.0,0.0) # BM3R1
cmap['ECK120029600'] = (0.0,0.0,0.0) # SrpR
cmap['ECK120033737'] = (0.0,0.0,0.0) # PhlF
cmap['L3S2P21'] = (0.0,0.0,0.0) # YFP

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
				region = [new_fwd_profile[new_start_bp:new_end_bp], 
				          new_rev_profile[new_start_bp:new_end_bp]]
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

#############################################################################
# PREDICTED profiles
#############################################################################

cir_logic = {}
cir_logic['AmtR'] = [0,1,1,1,0,1,1,1]
cir_logic['LitR'] = [0,0,1,1,1,1,1,1]
cir_logic['BM3R1'] = [0,0,0,0,1,1,1,1]
cir_logic['SrpR'] = [1,1,1,1,1,0,0,0]
cir_logic['PhlF'] = [1,1,0,0,0,1,1,1]
cir_logic['YFP'] = [0,0,1,1,1,0,0,0]

pred_exp = {}
pred_exp['AmtR'] = [0.149, 8.759, 15.869, 24.479, 0.149, 8.759, 15.869, 24.479]
pred_exp['LitR'] = [0.142, 0.142, 15.862, 15.862, 7.44, 7.44, 23.16, 23.16]
pred_exp['BM3R1'] = [0.013, 0.013, 0.013, 0.013, 7.311, 7.311, 7.311, 7.311]
pred_exp['SrpR'] = [26.14, 13.877, 16.538, 11.173, 8.532, 0.807, 1.199, 0.372]
pred_exp['PhlF'] = [5.282, 5.331, 0.364, 0.36, 0.334, 4.503, 4.98, 5.247]
pred_exp['YFP'] = [0.088, 0.093, 12.706, 12.706, 13.428, 0.104, 0.098, 0.098]

conv_facs = {}
conv_facs['AmtR'] = 12423.20455
conv_facs['LitR'] = 35289.74773
conv_facs['BM3R1'] = 618.8727273
conv_facs['SrpR'] = 5037.768182
conv_facs['PhlF'] = 6913.379545
conv_facs['YFP'] = 3050.0

pred_fpkms = {}
for r_key in pred_exp.keys():
	pred_fpkms[r_key] = np.array(pred_exp[r_key]) * conv_facs[r_key]

def load_frag_profile (filename):
	data = [[],[]]
	f = open(filename, 'rU')
	reader = csv.reader(f, delimiter='\t')
	for row in reader:
		if len(row) == 2:
			data[0].append(float(row[0]))
			data[1].append(float(row[1]))
	profile = np.zeros(1000)
	for idx in range(len(data[0])):
		if data[0][idx] < 1000:
			profile[idx] = data[1][idx]
	return profile/np.sum(profile)

frag_profiles = []
for s in range(8):
	frags = {}
	frag_profiles.append(frags)
	cur_frag_dist = load_frag_profile('./frag_profiles/frag_profile_'+str(s+1)+'.txt')
	frags['AmtR'] = fm.frag_factor_profile(989-149, frag_dist=cur_frag_dist)
	frags['LitR'] = fm.frag_factor_profile(2178-1398, frag_dist=cur_frag_dist)
	frags['BM3R1'] = fm.frag_factor_profile(3269-2513, frag_dist=cur_frag_dist)
	frags['SrpR'] = fm.frag_factor_profile(4247-3402, frag_dist=cur_frag_dist)
	frags['PhlF'] = fm.frag_factor_profile(5172-4396, frag_dist=cur_frag_dist)
	frags['YFP'] = fm.frag_factor_profile(6124-5242, frag_dist=cur_frag_dist)

pred_traces_on = []
pred_traces_off = []
for s in range(8):
	pred_trace_on = np.zeros(7124)
	pred_trace_off = np.zeros(7124)
	if cir_logic['AmtR'][s] == 1:
		pred_trace_on[149:989] = np.ones(989-149)*pred_fpkms['AmtR'][s]*frag_profiles[s]['AmtR']
	else:
		pred_trace_off[149:989] = np.ones(989-149)*pred_fpkms['AmtR'][s]*frag_profiles[s]['AmtR']
	if cir_logic['LitR'][s] == 1:
		pred_trace_on[1398:2178] = np.ones(2178-1398)*pred_fpkms['LitR'][s]*frag_profiles[s]['LitR']
	else:
		pred_trace_off[1398:2178] = np.ones(2178-1398)*pred_fpkms['LitR'][s]*frag_profiles[s]['LitR']
	if cir_logic['BM3R1'][s] == 1:
		pred_trace_on[2513:3269] = np.ones(3269-2513)*pred_fpkms['BM3R1'][s]*frag_profiles[s]['BM3R1']
	else:
		pred_trace_off[2513:3269] = np.ones(3269-2513)*pred_fpkms['BM3R1'][s]*frag_profiles[s]['BM3R1']
	if cir_logic['SrpR'][s] == 1:
		pred_trace_on[3402:4247] = np.ones(4247-3402)*pred_fpkms['SrpR'][s]*frag_profiles[s]['SrpR']
	else:
		pred_trace_off[3402:4247] = np.ones(4247-3402)*pred_fpkms['SrpR'][s]*frag_profiles[s]['SrpR']
	if cir_logic['PhlF'][s] == 1:
		pred_trace_on[4396:5172] = np.ones(5172-4396)*pred_fpkms['PhlF'][s]*frag_profiles[s]['PhlF']
	else:
		pred_trace_off[4396:5172] = np.ones(5172-4396)*pred_fpkms['PhlF'][s]*frag_profiles[s]['PhlF']
	if cir_logic['YFP'][s] == 1:
		pred_trace_on[5242:6124] = np.ones(6124-5242)*pred_fpkms['YFP'][s]*frag_profiles[s]['YFP']
	else:
		pred_trace_off[5242:6124] = np.ones(6124-5242)*pred_fpkms['YFP'][s]*frag_profiles[s]['YFP']
	pred_traces_on.append(pred_trace_on)
	pred_traces_off.append(pred_trace_off)

#############################################################################
# PLOT profiles
#############################################################################

def plot_act_vs_pred_traces (act_traces, pred_traces_on, pred_traces_off, dna_design, x_range, out_filename):
	fig = plt.figure(figsize=(5.0,3.6)) # 5,5 for SI 5.8,5.0
	gs = gridspec.GridSpec(9, 1, height_ratios=[1,1,1,1,1,1,1,1,0.8])
	# Plot the DNA
	dr = dpl.DNARenderer(scale=2, linewidth=0.75)
	ax = plt.subplot(gs[8])
	start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
	ax.set_xlim([0, 6124])
	ax.set_ylim([-10,8])
	ax.plot([0,6124], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')
	# Plot the traces
	for s in range(8):
		ax = plt.subplot(gs[s])
		act_trace = act_traces[s]
		ax.fill_between(range(x_range[1]), act_trace, np.zeros(len(act_trace)), facecolor=col_sense_rd, linewidth=0)
		ax.plot(range(len(pred_traces_on[s])), pred_traces_on[s], color=col_on_state, linewidth=1.0)
		ax.plot(range(len(pred_traces_on[s])), pred_traces_off[s], color=col_off_state, linewidth=1.0)
		ax.set_xlim([x_range[0],x_range[1]])
		ax.set_xticks([])
		ax.set_yscale('symlog', linthreshx=10000)
		ax.set_ylim([5, 8.0*10e5])
		ax.set_yticks([10e1, 10e3, 10e5])
		ax.tick_params(axis='y', which='major', labelsize=8, pad=1, length=2, width=0.5)
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(0.8)
	plt.subplots_adjust(hspace=.00, wspace=.00, left=.07, right=.99, top=.99, bottom=.01)
	fig.savefig(out_filename, transparent=True)
	plt.close('all')

def gff_to_dnaplotlib (gff, chrom):
	design = []
	for part_name in gff[chrom].keys():
		part_data = gff[chrom][part_name]
		p_name = part_name
		p_type = None
		p_opts = {}
		if part_data[0] == 'gene':
			p_type = 'CDS'
			p_opts['label'] = part_name
			p_opts['label_y_offset'] = -6.5
			p_opts['label_size'] = 8
			p_opts['label_style'] = 'italic'
		if part_data[0] == 'promoter':
			p_type = 'Promoter'
		if part_data[0] == 'terminator':
			p_type = 'Terminator'
		if p_type != None:
			p_start_bp = int(part_data[2])
			p_end_bp = int(part_data[3])
			if p_name in cmap.keys():
				p_opts['color'] = cmap[p_name]
			new_part = {'type':p_type,
				        'name':p_name,
				        'start':p_start_bp,
				        'end':p_end_bp,
				        'fwd':True,
				        'opts':p_opts}
			design.append(new_part)
	return sorted(design, key=lambda k: k['start']) 

def plot_tube_vs_flask_traces (tube_traces, flask_traces, dna_design, x_range, out_filename, text_size=8):
	fig = plt.figure(figsize=(5.0,3.6)) # 5,5 for SI 5.8,5.0
	gs = gridspec.GridSpec(9, 1, height_ratios=[1,1,1,1,1,1,1,1,0.8])
	# Plot the DNA
	dr = dpl.DNARenderer(scale=2, linewidth=0.75)
	ax = plt.subplot(gs[8])
	start, end = dr.renderDNA(ax, dna_design, dr.trace_part_renderers())
	ax.set_xlim([0, 6124])
	ax.set_ylim([-10,8])
	ax.plot([0,6124], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.axis('off')
	# Plot the traces
	for s in range(8):
		ax = plt.subplot(gs[s])
		flask_trace = flask_traces[s]
		tube_trace = tube_traces[s]
		ax.fill_between(range(x_range[1]), flask_trace, np.zeros(len(flask_trace)), facecolor=col_sense_rd, linewidth=0)
		ax.plot(range(x_range[1]), tube_trace, color=(0,0,0), linewidth=1, linestyle='-')
		ax.set_xlim([x_range[0],x_range[1]])
		ax.set_xticks([])
		ax.set_yscale('symlog', linthreshx=10000)
		ax.set_ylim([5, 10e5])
		ax.set_yticks([10e1, 10e3])
		ax.tick_params(axis='y', which='major', labelsize=text_size, pad=1, length=2, width=0.5)
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(0.8)
	plt.subplots_adjust(hspace=.00, wspace=.00, left=.07, right=.99, top=.99, bottom=.01)
	fig.savefig(out_filename, transparent=True)
	plt.close('all')

# Load design information and transcription profiles
gff = load_gff(DATA_PREFIX+'/gff/0x58v50.gff')
dna_design = gff_to_dnaplotlib(gff, '0x58v50')
norm_factors = load_norm_factors(RESULTS_PREFIX+'/norm.factors.matrix.txt')
tube_profiles = {}
for s in tube_samples:
	tube_profiles[s] = load_profiles(RESULTS_PREFIX+'/'+s+'/'+s+'.fwd.profiles.txt', RESULTS_PREFIX+'/'+s+'/'+s+'.rev.profiles.txt')
flask_profiles = {}
for s in flask_samples:
	flask_profiles[s] = load_profiles(RESULTS_PREFIX+'/'+s+'/'+s+'.fwd.profiles.txt', RESULTS_PREFIX+'/'+s+'/'+s+'.rev.profiles.txt')

# Calculate a conversion factor from FPKM to reads use tube_4 as high number of reads and YFP FPKM
YFP_reads = np.array(extract_profile_region(tube_profiles['tube_4'], '0x58v50', gff['0x58v50']['YFP'][2]+100, gff['0x58v50']['YFP'][3]-100)[0])/norm_factors['tube_4']
YFP_FPKM = 38977.39525
reads_FPKM_fac = YFP_FPKM/np.median(YFP_reads)

# Load the trace data
tube_traces = []
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_1'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_1']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_2'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_2']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_3'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_3']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_4'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_4']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_5'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_5']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_6'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_6']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_7'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_7']) * reads_FPKM_fac )
tube_traces.append( (np.array(extract_profile_region(tube_profiles['tube_8'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['tube_8']) * reads_FPKM_fac )
flask_traces = []
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_1'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_1']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_2'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_2']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_3'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_3']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_4'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_4']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_5'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_5']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_6'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_6']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_7'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_7']) * reads_FPKM_fac )
flask_traces.append( (np.array(extract_profile_region(flask_profiles['flask_8'], '0x58v50', 0, gff['0x58v50']['L3S2P21'][3])[0]) / norm_factors['flask_8']) * reads_FPKM_fac )
x_range = [0, gff['0x58v50']['L3S2P21'][3]]

# Plot the two figures
plot_act_vs_pred_traces(tube_traces, pred_traces_on, pred_traces_off, dna_design, x_range, OUTPUT_PREFIX+'/circuit_profile_tube.pdf')


#Make text bigger
for part in dna_design:
	if part['type'] == 'CDS':
		part['opts']['label_y_offset'] = -7.3
		part['opts']['label_size'] = 9.5

plot_tube_vs_flask_traces(tube_traces, flask_traces, dna_design, x_range, OUTPUT_PREFIX+'/circuit_profile_flask_vs_tube.pdf', text_size=9.5)

