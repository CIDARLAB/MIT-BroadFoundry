#!/usr/bin/env python
"""
	Workflow for the circuit data set - PART I
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
FASTQ_PREFIX = '/btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/FastqsConcatenated_Technical_Replicates/'
HOME_DIR = '/home/unix/tgorocho/genetic_debugger/'
DATA_PREFIX = HOME_DIR+'runs/circuit/data/'
RESULTS_PREFIX = '/broad/hptmp/tgorocho/circuit/'

# Split the designs into a set number of jobs (4 processors will be used)
print 'workflow_circuit_1.py INFO: Processing designs in these blocks:'
blocked_designs = [['1_1','2_1','3_1','4_1'],
                   ['5_1','6_1','7_1','8_1'], 
                   ['1_2','2_2','3_2','4_2'],
                   ['5_2','6_2','7_2','8_2']]

for block in blocked_designs:
	print 'Block', block

###############################################################################
# RUN THE WORKFLOW
###############################################################################

def create_dir_if_needed (directory):
	"""Create a directory if it doesn't exist
	"""
	if not os.path.exists(directory):
		try:
			os.makedirs(directory)
		except IOError as e:
			print("WARNING: Directory already exists.")

# Make sure the logs directory exists
create_dir_if_needed(RESULTS_PREFIX+'logs')

# For each design run eXpress (in parallel)
for designs in blocked_designs:
	err_log_name = 'circuit_1-'+'_'.join(designs)+'.err'
	out_log_name = 'circuit_1-'+'_'.join(designs)+'.out'
	cmd_to_run   = 'qsub -N cir1_'+designs[0]+'-'+designs[-1] + \
	               ' -e '+RESULTS_PREFIX+'logs/'+err_log_name + \
	               ' -o '+RESULTS_PREFIX+'logs/'+out_log_name + \
	               ' -b y -q gaag' + \
	               ' "python '+HOME_DIR+'bin/01_map_expression.py'+ \
	               ' -designs '+','.join(designs) + \
	               ' -fastq_prefix '+FASTQ_PREFIX + \
	               ' -data_prefix '+DATA_PREFIX + \
	               ' -results_prefix '+RESULTS_PREFIX + \
	               ' -clean_up Y"'
	print 'workflow_circuit_1.py RUNNING:', cmd_express_run
	subprocess.call(cmd_to_run, shell=True)
