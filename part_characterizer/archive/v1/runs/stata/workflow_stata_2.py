#!/usr/bin/env python
"""
	Workflow for the Stata data set - PART II
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
LIBRARY_FILE = DATA_PREFIX+'stata_library.txt'

# Which designs to create estimates for (85 to process all designs)
designs_to_process = [str(x) for x in range(1,85)]
# Remove designs with sequencing errors
temp_designs = []
for idx in range(len(designs_to_process)):
	if designs_to_process[idx] not in ['17', '18', '33', '45', '57', '69', '76', '81']:
		temp_designs.append(designs_to_process[idx])
designs_to_process = temp_designs

###############################################################################
# RUN THE WORKFLOW
###############################################################################

cmd_express_run = 'python '+HOME_DIR+'bin/02_normalize.py'+ \
				      ' -library '+LIBRARY_FILE + \
	                  ' -designs '+','.join(designs) + \
	                  ' -results_prefix '+RESULTS_PREFIX + \
	                  ' -normaliser edgeR'
print 'workflow_stata_2.py RUNNING:', cmd_express_run
subprocess.call(cmd_express_run, shell=True)

cmd_express_run = 'python '+HOME_DIR+'bin/03_part_performance.py'+ \
				      ' -library '+LIBRARY_FILE + \
	                  ' -designs '+','.join(designs) + \
	                  ' -measure fpkm' + \
	                  ' -data_prefix '+RESULTS_PREFIX + \
	                  ' -results_prefix '+RESULTS_PREFIX+'part_performance/'
print 'workflow_stata_2.py RUNNING:', cmd_express_run
subprocess.call(cmd_express_run, shell=True)
