#!/usr/bin/env python
"""
Generate job shell scripts for GridEngine cluster
=================================================
	Create the shell scripts to run each pool spearately on the
	Broad GridEngine cluster.
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
run_prefix = home_prefix + 'runs/2015-09_nif/'
fastq_prefix = '/btl/projects/Foundry/Yongjin/Tom/nif_fastq/'
regex_prefix   = run_prefix + 'data/regex_refs/'
results_prefix = run_prefix + 'results/'
bin_path       = home_prefix + 'bin/'

# The indexes to use for the extraction barcodes (normally 1 and 2, unless other tracking barcodes included)
R1_fwd_bc_idx = '1'
R1_rev_bc_idx = '2'
R2_fwd_bc_idx = '-2'
R2_rev_bc_idx = '-1'
other_bc_idxs = '0'

# The names of the raw read files and the length of the forward and reverse primer.
run_data = {}
run_data['nifB'] = ['nifB_R1.fastq', 'nifB_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
run_data['nifD'] = ['nifD_R1.fastq', 'nifD_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifE'] = ['nifE_R1.fastq', 'nifE_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifF'] = ['nifF_R1.fastq', 'nifF_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifH'] = ['nifH_R1.fastq', 'nifH_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifJ'] = ['nifJ_R1.fastq', 'nifJ_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
run_data['nifK'] = ['nifK_R1.fastq', 'nifK_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifM'] = ['nifM_R1.fastq', 'nifM_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifN'] = ['nifN_R1.fastq', 'nifN_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifQ'] = ['nifQ_R1.fastq', 'nifQ_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifS'] = ['nifS_R1.fastq', 'nifS_R2.fastq', '20', '20', R2_fwd_bc_idx, R2_rev_bc_idx, other_bc_idxs]
run_data['nifU'] = ['nifU_R1.fastq', 'nifU_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
run_data['nifV'] = ['nifV_R1.fastq', 'nifV_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
run_data['nifW'] = ['nifW_R1.fastq', 'nifW_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
run_data['nifY'] = ['nifY_R1.fastq', 'nifY_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
run_data['nifZ'] = ['nifZ_R1.fastq', 'nifZ_R2.fastq', '20', '20', R1_fwd_bc_idx, R1_rev_bc_idx, other_bc_idxs]
# These nifW pools we treat differently (read starts with barcode)
run_data['nifW_old'] = ['nifW_old_R1.fastq', 'nifW_old_R2.fastq', '0', '0', '1', '-1', other_bc_idxs]
run_data['nifW_new'] = ['nifW_new_R1.fastq', 'nifW_new_R2.fastq', '0', '0', '1', '-1', other_bc_idxs]

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
