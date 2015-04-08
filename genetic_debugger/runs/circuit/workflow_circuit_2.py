#!/usr/bin/env python
"""
	Run eXpress expression estimates workflow for the circuit data set
"""
#	Copyright (C) 2014 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import os
import subprocess

###############################################################################
# PARAMETERS FOR THE WORKFLOW
###############################################################################

# Paths to data and outputs
HOME_DIR = '/home/unix/tgorocho/genetic_debugger/'
DATA_PREFIX = HOME_DIR+'runs/circuit/data/'
RESULTS_PREFIX = '/broad/hptmp/tgorocho/circuit/'

designs_to_process = ['1_1', '2_1', '3_1', '4_1', '5_1', '6_1', '7_1', '8_1', 
                      '1_2', '2_2', '3_2', '4_2', '5_2', '6_2', '7_2', '8_2']

###############################################################################
# RUN THE WORKFLOW
###############################################################################

run_express_normalizer(designs_to_process, 
	                   RESULTS_PREFIX + 'gene_exp/',
	                   RESULTS_PREFIX + 'transcript_exp/',
	                   normaliser='DESeq')

run_express_normalizer(designs_to_process, 
	                   RESULTS_PREFIX + 'gene_exp/',
	                   RESULTS_PREFIX + 'transcript_exp/',
	                   normaliser='edgeR')
