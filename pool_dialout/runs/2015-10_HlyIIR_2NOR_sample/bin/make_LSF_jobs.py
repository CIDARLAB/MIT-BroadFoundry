#!/usr/bin/env python
"""
Generate job shell scripts for LSF cluster
=================================================
	Create the shell scripts to run each pool spearately on the
	Broad LSF cluster.
"""
from __future__ import print_function, division
import os
import sys
import string

__author__  = 'Thomas E. Gorochowski, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

# Prefix for saving shell scripts to
save_to = '../'

# Prefixes for files
home_prefix = '/home/unix/tgorocho/pool_dialout/'
run_prefix = home_prefix + 'runs/2015-10_HlyIIR_2NOR_sample/'
fastq_prefix = '/btl/projects/Foundry/Lauren/barcodes_Sept2015/'
regex_prefix   = run_prefix + 'data/regex_refs/'
results_prefix = run_prefix + 'results/'
bin_path       = home_prefix + 'bin/'

# The indexes to use for the extraction barcodes (normally 1 and 2, unless other tracking barcodes included)
R1_fwd_bc_idx = '2'
R1_rev_bc_idx = '3'
R2_fwd_bc_idx = ''
R2_rev_bc_idx = ''
other_bc_idxs = '1'

# The names of the raw read files and the length of the forward and reverse primer.
run_data = {}
run_data['HlyIIR_2NOR'] = ['spl1_S1_L001_R1_001_sample.fastq', 'spl1_S1_L001_R2_001_sample.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]

# Build the part analysis script iteratively
part_analysis_cmd = ''

# Create the script for each pool
for name in run_data.keys():
	# Standard places to put things
	stdout_file = results_prefix + name + '_stdout.txt'
	stderr_file = results_prefix + name + '_stderr.txt'
	bin_file = bin_path + 'perfect_dialout.py'
	regex_file = regex_prefix + name + '_regexs.txt'
	# Build the command to use GridEngine
	cmd =  ''
	cmd += 'bsub -e ' + stderr_file
	cmd += ' -o ' + stdout_file
	cmd += ' -q forest -N '
	cmd += ' "python ' + bin_file
	cmd += ' ' + regex_file
	cmd += ' ' + fastq_prefix + run_data[name][0]
	cmd += ' ' + fastq_prefix + run_data[name][1]
	cmd += ' ' + run_data[name][2]
	cmd += ' ' + run_data[name][3]
	cmd += ' ' + run_data[name][4]
	cmd += ' ' + run_data[name][5]
	cmd += ' ' + run_data[name][6]
	cmd += ' ' + results_prefix + name + '_"'
	# Save the script file
	f_out = open(save_to + 'run_'+name+'.sh', 'w')
	f_out.write(cmd)
	f_out.close()
	# Add part analysis commmand
	part_analysis_cmd += 'python ../../bin/dialout_part_analysis.py ./results/' + name + '_dialout_designs.csv  ./results/' + name + '_dialout_part_analysis.csv\n'

# Script to run everything
f_out = open(save_to + 'run_all_part_analysis.sh', 'w')
f_out.write(part_analysis_cmd)
f_out.close()

# Script to run everything
f_out = open(save_to + 'run_all.sh', 'w')
for name in sorted(run_data.keys()):
	f_out.write('sh  run_'+name+'.sh\n')
f_out.close()
