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
run_prefix = home_prefix + 'runs/2015-03-01_circuits/'
fastq_prefix = '/btl/data/MiSeq0/runs/lauren/150227_M03102_0065_000000000-AE6TR/Data/Intensities/BaseCalls/'
regex_prefix   = run_prefix + 'data/regex_refs/'
results_prefix = run_prefix + 'results/'
bin_path       = home_prefix + 'bin/'

# The indexes to use for the extraction barcodes (normally 1 and 2, unless other tracking barcodes included)
fwd_barcode_idx = 2
rev_barcode_idx = 3

run_data = {}
run_data['AmeR_NOT'] = ['AmeR-NOT-tagged-57_S1_L001_R1_001.fastq',      'AmeR-NOT-tagged-57_S1_L001_R2_001.fastq', '20', '20']
run_data['AmtR_NOT'] = ['AmtR-NOT-tagged-726_S2_L001_R1_001.fastq',     'AmtR-NOT-tagged-726_S2_L001_R2_001.fastq', '20', '20']
run_data['BetI_NOT'] = ['BetI-NOT-tagged-630_S3_L001_R1_001.fastq',     'BetI-NOT-tagged-630_S3_L001_R2_001.fastq', '20', '20']
run_data['BM3R1_NOT'] = ['BM3R1-NOT-tagged-954_S4_L001_R1_001.fastq',   'BM3R1-NOT-tagged-954_S4_L001_R2_001.fastq', '20', '20']
run_data['HlyIIR_NOT'] = ['HlyIIR-NOT-tagged-190_S5_L001_R1_001.fastq', 'HlyIIR-NOT-tagged-190_S5_L001_R2_001.fastq', '20', '20']
run_data['IcaRA_NOT'] = ['IcaRA-NOT-tagged-332_S6_L001_R1_001.fastq',   'IcaRA-NOT-tagged-332_S6_L001_R2_001.fastq', '20', '20']
run_data['LitR_NOT'] = ['LitR-NOT-tagged-393_S7_L001_R1_001.fastq',     'LitR-NOT-tagged-393_S7_L001_R2_001.fastq', '20', '20']
run_data['PhlF_NOT'] = ['PhlF-NOT-tagged-52_S8_L001_R1_001.fastq',      'PhlF-NOT-tagged-52_S8_L001_R2_001.fastq', '20', '20']
run_data['QacR_NOT'] = ['QacR-NOT-tagged-100_S9_L001_R1_001.fastq',     'QacR-NOT-tagged-100_S9_L001_R2_001.fastq', '20', '20']
run_data['SrpR_NOT'] = ['SrpR-NOT-tagged-375_S10_L001_R1_001.fastq',    'SrpR-NOT-tagged-375_S10_L001_R2_001.fastq', '20', '20']

run_data['AmeR_2NOR'] = ['AmeR-2NOR-tagged-741_S11_L001_R1_001.fastq',     'AmeR-2NOR-tagged-741_S11_L001_R2_001.fastq', '20', '20']
run_data['AmtR_2NOR'] = ['AmtR-2NOR-tagged-908_S12_L001_R1_001.fastq',     'AmtR-2NOR-tagged-908_S12_L001_R2_001.fastq', '20', '20']
run_data['BetI_2NOR'] = ['BetI-2NOR-tagged-34_S13_L001_R1_001.fastq',      'BetI-2NOR-tagged-34_S13_L001_R2_001.fastq', '20', '20']
run_data['BM3R1_2NOR'] = ['BM3R1-2NOR-tagged-236_S14_L001_R1_001.fastq',   'BM3R1-2NOR-tagged-236_S14_L001_R2_001.fastq', '20', '20']
run_data['HlyIIR_2NOR'] = ['HlyIIR-2NOR-tagged-869_S15_L001_R1_001.fastq', 'HlyIIR-2NOR-tagged-869_S15_L001_R2_001.fastq', '20', '20']
run_data['IcaRA_2NOR'] = ['IcaRA-2NOR-tagged-960_S16_L001_R1_001.fastq',   'IcaRA-2NOR-tagged-960_S16_L001_R2_001.fastq', '20', '20']
run_data['LitR_2NOR'] = ['LitR-2NOR-tagged-930_S17_L001_R1_001.fastq',     'LitR-2NOR-tagged-930_S17_L001_R2_001.fastq', '20', '20']
run_data['PhlF_2NOR'] = ['PhlF-2NOR-tagged-214_S18_L001_R1_001.fastq',     'PhlF-2NOR-tagged-214_S18_L001_R2_001.fastq', '20', '20']
run_data['QacR_2NOR'] = ['QacR-2NOR-tagged-367_S19_L001_R1_001.fastq',     'QacR-2NOR-tagged-367_S19_L001_R2_001.fastq', '20', '20']
run_data['SrpR_2NOR'] = ['SrpR-2NOR-tagged-426_S20_L001_R1_001.fastq',     'SrpR-2NOR-tagged-426_S20_L001_R2_001.fastq', '20', '20']


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
