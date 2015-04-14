#!/usr/bin/env python
"""
	Workflow for the circuit data set - PART II
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
LIBRARY_FILE = DATA_PREFIX+'circuit_library.txt'

designs_to_process = ['1_1', '2_1', '3_1', '4_1', '5_1', '6_1', '7_1', '8_1', 
                      '1_2', '2_2', '3_2', '4_2', '5_2', '6_2', '7_2', '8_2']

###############################################################################
# RUN THE WORKFLOW
###############################################################################

cmd_express_run = 'python '+HOME_DIR+'bin/02_normalize.py'+ \
				      ' -library '+LIBRARY_FILE + \
	                  ' -designs '+','.join(designs) + \
	                  ' -results_prefix '+RESULTS_PREFIX + \
	                  ' -normaliser edgeR'
print 'workflow_circuit_2.py RUNNING:', cmd_express_run
subprocess.call(cmd_express_run, shell=True)

cmd_express_run = 'python '+HOME_DIR+'bin/03_part_performance.py'+ \
				      ' -library '+LIBRARY_FILE + \
	                  ' -designs '+','.join(designs) + \
	                  ' -measure fpkm' + \
	                  ' -data_prefix '+RESULTS_PREFIX + \
	                  ' -results_prefix '+RESULTS_PREFIX+'part_performance/'
print 'workflow_circuit_2.py RUNNING:', cmd_express_run
subprocess.call(cmd_express_run, shell=True)
