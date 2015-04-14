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
designs_to_process = ['7', '21', '27', '36', '42', '52', '60', '68', '74']
print 'workflow_stata_express.py INFO: Processing designs:', designs_to_process

blocked_designs = [['7', '21'], ['27', '36'], ['42', '52'], ['60', '68'], ['74']]

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
