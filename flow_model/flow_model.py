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

    Note: this is a deterministic model, not taking stochastic aspects
    into consideration. As the data we generally compare to (RNA-seq)
    is averaged across the entire populations, stochastic effects will
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

def generate_site_model (gcl, variant, part_start_idx, part_end_idx, site_len=10):
	# Sites along the DNA molecule (0 = +ve strand, 1 = -ve strand)
	# Sites can have the following classifications:
	# 	0 - non-coding DNA
	# 	1 - coding DNA (CDS)
	# 	2 - promoter 
	# 	3 - terminator
	sites = [[],[]]
	# Abstract the DNA sequence from the library into a site model using types above



	# Rates for transitions between sites
	rate_fwd = [[],[]]
	rate_rev = [[],[]]
	rate_int = [[],[]]
	rate_ext = [[],[]]
	rates = [rate_fwd, rate_rev, rate_int, rate_ext]
	return (sites, rates)

def flow_derivative(y, time, rate_fwd, rate_rev, rate_int, rate_ext):
	# vector to hold output
	out = np.zeros(y)
	# Output vetor is actually a concatentation of [+ve,-ve] strand data
	s_num = np.size(out, 0)/2.0
	# Cycle through each pair (+ve,-ve) sites and calculate derivative
	for i in range(s_num):
		if i == 0:
			# Handle start points
			out[i] = y[i]*rate_int[i] + y[i+1]*rate_rev[i+1] - y[i]*rate_fwd[i] - y[i]*rate_rev[i] - y[i]*rate_ext[i]
		elif i == s_num-1:
			# Handle end points
			out[i] = y[i]*rate_int[i] + y[i-1]*rate_fwd[i-1] - y[i]*rate_fwd[i] - y[i]*rate_rev[i] - y[i]*rate_ext[i]
		else:
			# Handle all other positions
			# y[i]
			y_i_f = y[i]
			y_i_r = y[i+s_num]
			y_i_all = y_i_f + y_i_r
			# y[i-1]
			y_i_prev_f = y[i-1]
			y_i_prev_r = y[i-1+s_num]
			y_i_prev_all = y_i_f + y_i_r
			# y[i+1]
			y_i_next_f = y[i-1]
			y_i_next_r = y[i-1+s_num]
			y_i_next_all = y_i_f + y_i_r
			# Update +ve strand
			out[i] = y_i_all*rate_int[i] + y_i_prev_all*rate_fwd[i-1] + y_i_next_all*rate_rev[i+1] - y_i_f*rate_fwd[i] - y_i_f*rate_rev[i] - y_i_f*rate_ext[i]
			# Update -ve strand
			out[i+s_num] = 0.0
	return out

def run_flow_model (model, t_vec):
	# TODO: add steady state checks (continue until stability reached)
	sites = model[0]
	rates = model[1]
	init_cond = np.zeros(len(sites)*2)
	y, info = odeint(flow_derivative, init_cond, time_vec, 
		             args=(rates[0], rates[1], rates[2], rates[3]), 
		             full_output=True)
	return y, info


# Test the model
model = generate_site_model(nifs, '1', 1, -1, site_len=10)
t_vec = np.linspace(0, 100, 40)
y, info = run_flow_model(model, t_vec)


