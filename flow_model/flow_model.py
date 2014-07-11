#!/usr/bin/env python
"""
Flow Model
==========

    Implementation of a flow-based model of gene transcription similar
    to that of T. Tuller's work:

    Edri et al. The RNA Polymerase Flow Model of Gene Transcription. 
    IEEE Transactions on Biomedical Circuits and Systems, 8:1, 2014.

    This enables the assessment of RNAP movement and final steady state 
    densities along a stretch of DNA taking into consideration interactions 
    between RNAPs moving potentially in different directions along the 
    molecule. Input is a region of a GeneClusterLibrary variant and output
    is the steady-state RNAP density distribution along the DNA.

    Note: this is a deterministic model not taking stochastic aspects
    into consideration. As the data we generally compare to (RNA-seq)
    is averaged across entire populations, stochastic effects will
    become smoothed (although may have significant effect at the single
    cell level).
"""
#    Flow Model
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
from scipy.integrate import odeint

FLOW_MODEL_PARAMS = {'fwd_rate' : 0.2,
                     'rev_rate' : 0.01,
                     'int_rate' : 0.4,
                     'ext_rate' : 0.9}

def generate_site_model (gcl, variant, part_start_idx, part_end_idx, site_len=20):
	# Abstract the DNA sequence from the library into a site model using types below
	start_bp = 0
	end_bp = 1
	# Generate the sites (initially classify as generic DNA)
	num_of_sites = int((end_bp-start_bp)/site_len)
	if (end_bp-start_bp) % site_len != 0:
		num_of_sites += 1
	# Sites along the DNA molecule (0 = +ve strand, 1 = -ve strand)
	# Sites can have the following classifications:
	# 	0 - generic DNA
	# 	1 - promoter
	# 	2 - coding DNA (CDS)
	# 	3 - terminator
	sites = [[0]*num_of_sites, [0]*num_of_sites]
	# Cycle through all parts 
	for cur_idx in range(part_start_idx, part_end_idx+1):
		# Extract the part and find site it belongs
		cur_el = gcl.variants[variant]['part_list'][cur_idx]
		cur_type = gcl.parts[cur_el['part_name']]['type']
		cur_start_bp = cur_el['seq_idx']
		cur_end_bp = cur_el['seq_idx'] + cur_el['seq_len']
		if cur_el['dir'] == 'R':
			cur_start_bp = temp
			cur_start_bp = cur_end_bp
			cur_end_bp = temp
		site_start = int((cur_start_bp-start_bp)/site_len)
		site_end = int((cur_end_bp-start_bp)/site_len)
		# Check part type and how to handle
		if cur_type == 'Promoter':
			# It's a promoter - use end bp
			if cur_el['dir'] == 'F':
				sites[0][site_end] = 1
			else:
				sites[1][site_end] = 1
		if cur_type == 'Terminator':
			# It's a terminator - use start bp
			if cur_el['dir'] == 'F':
				sites[0][site_start] = 3
			else:
				sites[1][site_start] = 3
		if cur_type == 'CDS':
			# It's a coding region - extend across (make sure not to overwrite promoter)
			if cur_el['dir'] == 'F':
				for s in range(site_start, site_end+1):
					if sites[0][s] == 0:
						sites[0][s] = 2
			else:
				for s in range(site_end, site_start+1):
					if sites[1][s] == 0:
						sites[1][s] = 2
	# Rates for transitions between sites (use defaults for fwd and rev rates - maybe change later)
	rate_fwd = [[FLOW_MODEL_PARAMS['fwd_rate']]*num_of_sites,[FLOW_MODEL_PARAMS['fwd_rate']]*num_of_sites]
	rate_rev = [[FLOW_MODEL_PARAMS['rev_rate']]*num_of_sites,[FLOW_MODEL_PARAMS['rev_rate']]*num_of_sites]
	rate_int = [[0.0]*num_of_sites,[0.0]*num_of_sites]
	rate_ext = [[0.0]*num_of_sites,[0.0]*num_of_sites]
	# Update the initation and termination rates based on site structure





	# TODO in progress





	rates = [rate_fwd, rate_rev, rate_int, rate_ext]
	return (sites, rates)

def flow_derivative(y, time, rate_fwd, rate_rev, rate_int, rate_ext):
	# vector to hold output
	out = np.zeros(y)
	# Output vetor is actually a concatentation of [+ve,-ve] strand data
	s_num = int(np.size(out, 0)/2)
	# Cycle through each pair (+ve,-ve) sites and calculate derivative
	for i in range(s_num):
		# y[i]
		y_i_f = y[i]
		y_i_r = y[i+s_num]
		y_i_all = y_i_f + y_i_r
		if i == 0:
			# Handle start points
			# y[next]
			y_i_next_f = y[i+1]
			y_i_next_r = y[i+1+s_num]
			y_i_next_all = y_i_f + y_i_r
			# Update +ve strand
			out[i] = (  rate_int[i]*(1.0-y_i_all)
                      + y_i_next_f*rate_rev[i+1]*(1.0-y_i_all)
                      - y_i_f*rate_fwd[i]*(1.0-y_i_next_all)
                      - y_i_f*rate_rev[i]
                      - rate_ext[i]*y_i_f )
			# Update -ve strand
			out[i+s_num] = (  rate_int[i+s_num]*(1.0-y_i_all)
                      + y_i_next_r*rate_fwd[i+1+s_num]*(1.0-y_i_all)
                      - y_i_r*rate_fwd[i+s_num]
                      - y_i_r*rate_rev[i+s_num]*(1.0-y_i_next_all)
                      - rate_ext[i+s_num]*y_i_r )
		elif i == s_num-1:
			# Handle end points
			# y[prev]
			y_i_prev_f = y[i-1]
			y_i_prev_r = y[i-1+s_num]
			y_i_prev_all = y_i_f + y_i_r
			# Update +ve strand
			out[i] = (  rate_int[i]*(1.0-y_i_all)
                      + y_i_prev_f*rate_fwd[i-1]*(1.0-y_i_all)
                      - y_i_f*rate_fwd[i]
                      - y_i_f*rate_rev[i]*(1.0-y_i_prev_all)
                      - rate_ext[i]*y_i_f )
			# Update -ve strand
			out[i+s_num] = (  rate_int[i+s_num]*(1.0-y_i_all)
                      + y_i_prev_r*rate_rev[i-1+s_num]*(1.0-y_i_all)
                      - y_i_r*rate_fwd[i+s_num]*(1.0-y_i_prev_all)
                      - y_i_r*rate_rev[i+s_num]
                      - rate_ext[i+s_num]*y_i_r )
		else:
			# Handle all other positions
			# y[prev]
			y_i_prev_f = y[i-1]
			y_i_prev_r = y[i-1+s_num]
			y_i_prev_all = y_i_f + y_i_r
			# y[next]
			y_i_next_f = y[i+1]
			y_i_next_r = y[i+1+s_num]
			y_i_next_all = y_i_f + y_i_r
			# Update +ve strand
			out[i] = (  rate_int[i]*(1.0-y_i_all)
                      + y_i_prev_f*rate_fwd[i-1]*(1.0-y_i_all)
                      + y_i_next_f*rate_rev[i+1]*(1.0-y_i_all)
                      - y_i_f*rate_fwd[i]*(1.0-y_i_next_all)
                      - y_i_f*rate_rev[i]*(1.0-y_i_prev_all)
                      - rate_ext[i]*y_i_f )
			# Update -ve strand
			out[i+s_num] = (  rate_int[i+s_num]*(1.0-y_i_all)
                      + y_i_prev_r*rate_rev[i-1+s_num]*(1.0-y_i_all)
                      + y_i_next_r*rate_fwd[i+1+s_num]*(1.0-y_i_all)
                      - y_i_r*rate_fwd[i+s_num]*(1.0-y_i_prev_all)
                      - y_i_r*rate_rev[i+s_num]*(1.0-y_i_next_all)
                      - rate_ext[i+s_num]*y_i_r )
	return out

def run_flow_model (site_model, t_vec):
	# TODO: add steady state checks (continue until stability reached)
	# Make model suitable for standard ODE solver (flatten structure)
	sites = site_model[0][0] + site_model[0][1]
	rates[0] = site_model[1][0][0] + site_model[1][0][1] # FWD
	rates[1] = site_model[1][1][0] + site_model[1][1][1] # REV
	rates[2] = site_model[1][2][0] + site_model[1][2][1] # INT
	rates[3] = site_model[1][3][0] + site_model[1][3][1] # EXT
	# Set up solver and run
	init_cond = np.zeros(sites)
	y, info = odeint(flow_derivative, init_cond, time_vec, 
		             args=(rates[0], rates[1], rates[2], rates[3]), 
		             full_output=True)
	# Return the trajectory
	return y, info

# Test the model
model = generate_site_model(nifs, '1', 1, -1, site_len=10)
t_vec = np.linspace(0, 100, 40)
y, info = run_flow_model(model, t_vec)
