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
FASTQ_PREFIX = '/btl/projects/seq/foundry/voight_96/'
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

# Split the designs into a set number of jobs
print 'workflow_stata_express.py INFO: Processing designs in these blocks:'
designs_per_block = int(len(designs_to_process)/CONCURRENT_JOBS)
blocked_designs = []
for b_idx in range(CONCURRENT_JOBS):
	start_idx = b_idx*designs_per_block
	if b_idx < CONCURRENT_JOBS-1:
		blocked_designs.append(designs_to_process[start_idx:(start_idx+designs_per_block)])
		print 'Block', b_idx, '=', designs_to_process[start_idx:(start_idx+designs_per_block)]
	else:
		blocked_designs.append(designs_to_process[start_idx:])
		print 'Block', b_idx, '=', designs_to_process[start_idx:]

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
	print 'workflow_stata_express.py RUNNING:', cmd_express_run
	subprocess.call(cmd_express_run, shell=True)
