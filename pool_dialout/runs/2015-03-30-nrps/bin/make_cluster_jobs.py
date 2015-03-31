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
home_prefix = '/home/unix/tgorocho/pool_mapping/'
run_prefix = home_prefix + 'runs/2015-03-30-nrps/'
fastq_prefix = '/btl/data/MiSeq0/runs/Fang/150327_M03102_0078_000000000-AE72W/Data/Intensities/BaseCalls/'
regex_prefix   = run_prefix + 'data/regex_refs/'
results_prefix = run_prefix + 'results/'
bin_path       = home_prefix + 'bin/'

# The indexes to use for the extraction barcodes (normally 1 and 2, unless other tracking barcodes included)
fwd_barcode_idx = 1
rev_barcode_idx = 2

# The names of the raw read files and the length of the forward and reverse primer.
run_data = {}
run_data['BN'] = ['BN-tagged-332_S6_L001_R1_001.fastq',         'BN-tagged-332_S6_L001_R2_001.fastq', '18', '18']
run_data['BP'] = ['BP-tagged-954_S4_L001_R1_001.fastq',         'BP-tagged-954_S4_L001_R2_001.fastq', '18', '18']
run_data['CN'] = ['CN-tagged-393_S7_L001_R1_001.fastq',         'CN-tagged-393_S7_L001_R2_001.fastq', '18', '18']
run_data['FN'] = ['FN-tagged-52_S8_L001_R1_001.fastq',          'FN-tagged-52_S8_L001_R2_001.fastq', '18', '18']
run_data['FP'] = ['FP-tagged-190_S5_L001_R1_001.fastq',         'FP-tagged-190_S5_L001_R2_001.fastq', '18', '18']
run_data['Hypo'] = ['Hypo-tagged-57_S1_L001_R1_001.fastq',      'Hypo-tagged-57_S1_L001_R2_001.fastq', '18', '18']
run_data['Sfp'] = ['Sfp-tagged-630_S3_L001_R1_001.fastq',       'Sfp-tagged-630_S3_L001_R2_001.fastq', '18', '18']
run_data['Trans'] = ['Trans-tagged-726_S2_L001_R1_001.fastq',   'Trans-tagged-726_S2_L001_R2_001.fastq', '18', '18']


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
	cmd += 'qsub -N ' + name
	cmd += ' -e ' + stderr_file
	cmd += ' -o ' + stdout_file
	cmd += ' -b y -q gaag'
	cmd += ' "python ' + bin_file
	cmd += ' ' + regex_file
	cmd += ' ' + fastq_prefix + run_data[name][0]
	cmd += ' ' + fastq_prefix + run_data[name][1]
	cmd += ' ' + run_data[name][2]
	cmd += ' ' + run_data[name][3]
	cmd += ' ' + str(fwd_barcode_idx)
	cmd += ' ' + str(rev_barcode_idx)
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
