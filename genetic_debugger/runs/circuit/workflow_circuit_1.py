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
FASTQ_PREFIX = '/btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/FastqsConcatenated_Technical_Replicates/'
HOME_DIR = '/home/unix/tgorocho/genetic_debugger/'
DATA_PREFIX = HOME_DIR+'runs/circuit/data/'
RESULTS_PREFIX = '/broad/hptmp/tgorocho/circuit/'

# Split the designs into a set number of jobs (4 processors will be used)
print 'rsem_workflow.py INFO: Processing designs in these blocks:'
blocked_designs = [['1_1','2_1','3_1','4_1'],
                   ['5_1','6_1','7_1','8_1'], 
                   ['1_2','2_2','3_2','4_2'],
                   ['5_2','6_2','7_2','8_2']]

# For testing
#blocked_designs = [['1_1'],
#                   ['2_1']]
#designs_to_process = ['1_1', '2_1']

for block in blocked_designs:
	print 'Block', block

###############################################################################
# RUN THE WORKFLOW
###############################################################################

def create_dir_if_needed (directory):
	"""Create a directory if it doesn't exist
	"""
	if not os.path.exists(directory):
		os.makedirs(directory)

# Make sure the logs directory exists
create_dir_if_needed(RESULTS_PREFIX+'logs')

# For each design run eXpress (in parallel)
for designs in blocked_designs:
	err_log_name = 'express_run-'+'_'.join(designs)+'.err'
	out_log_name = 'express_run-'+'_'.join(designs)+'.out'
	cmd_express_run = 'qsub -N express_run_'+designs[0]+'-'+designs[-1] + \
	                  ' -e '+RESULTS_PREFIX+'logs/'+err_log_name + \
	                  ' -o '+RESULTS_PREFIX+'logs/'+out_log_name + \
	                  ' -b y -q gaag' + \
	                  ' "python '+HOME_DIR+'01_estimate_expression/express_run_bowtie2.py'+ \
	                  ' -designs '+','.join(designs) + \
	                  ' -fastq_prefix '+FASTQ_PREFIX + \
	                  ' -data_prefix '+DATA_PREFIX + \
	                  ' -results_prefix '+RESULTS_PREFIX + \
	                  ' -clean_up Y"'
	print 'workflow_circuit_express.py RUNNING:', cmd_express_run
	subprocess.call(cmd_express_run, shell=True)
