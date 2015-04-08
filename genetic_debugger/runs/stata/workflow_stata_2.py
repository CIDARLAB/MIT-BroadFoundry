#!/usr/bin/env python
"""
	Run eXpress expression estimates workflow for the Stata data set
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
DATA_PREFIX = HOME_DIR+'runs/stata/data/'
RESULTS_PREFIX = '/broad/hptmp/tgorocho/stata/'

# Number of concurrent jobs to run (WARNING: lots of disk activity required)
CONCURRENT_JOBS = 10

# Which designs to create estimates for (85 to process all designs)
designs_to_process = [str(x) for x in range(1,85)]
# Remove designs with sequencing errors
temp_designs = []
for idx in range(len(designs_to_process)):
	if designs_to_process[idx] not in ['17', '18', '33', '45', '57', '69', '76', '81']:
		temp_designs.append(designs_to_process[idx])
designs_to_process = temp_designs
print 'workflow_stata_express.py INFO: Processing designs:', designs_to_process

###############################################################################
# RUN THE WORKFLOW
###############################################################################

run_express_normalizer(designs_to_process,
	                   RESULTS_PREFIX + 'express_gene_exp/',
	                   RESULTS_PREFIX + 'express_exp/',
	                   normaliser='DESeq')
run_express_normalizer(designs_to_process,
	                   RESULTS_PREFIX + 'express_gene_exp/',
	                   RESULTS_PREFIX + 'express_exp/',
	                   normaliser='edgeR')

