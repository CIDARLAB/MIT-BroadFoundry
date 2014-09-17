#!/usr/bin/env python
"""
Analysis of the frgamentation distributions for the Stata nif clusters
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import fragmentation_model as fm
import matplotlib.pyplot as plt
import csv
import numpy as np
import shutil

def read_frag_dist (in_filename):
	reads_dict = {}
	# Load the file data
	file_reader = csv.reader(open(in_filename, 'rU'), delimiter='\t')
	at_data = False
	for row in file_reader:
		if at_data == False:
			if len(row) > 0 and row[0] == 'insert_size':
				at_data = True
		else:
			if len(row) == 2:
				reads_dict[int(row[0])] = float(row[1])
	reads = np.zeros(max(reads_dict.keys())+1)
	for i in range(len(reads+1)):
		if i in reads_dict.keys():
			reads[i] = reads_dict[i]
	reads_norm = np.array(reads)/np.sum(reads)
	return reads, reads_norm

def resave_frag_dists (convert_filename, input_path, output_path):
	# Load the conversion mapping to variant name
	convert_map = {}
	file_reader = csv.reader(open(convert_filename, 'rU'), delimiter=',')
	# Ignore the header
	header = next(file_reader)
	for row in file_reader:
		if len(row) == 2:
			convert_map[row[0]] = row[1]
	for f in convert_map.keys(): 
		cur_from_file = input_path + f + '_1.stats'
		cur_to_file = output_path + 'variant_' + convert_map[f] + '_1.stats'
		shutil.copy2(cur_from_file, cur_to_file)
		cur_from_file = input_path + f + '_2.stats'
		cur_to_file = output_path + 'variant_' + convert_map[f] + '_2.stats'
		shutil.copy2(cur_from_file, cur_to_file)

def generate_frag_data (frag_dist, mrna_lens, no_end=False):
	# Generate range of profiles for differing lengths mRNAs
	profiles = []
	peaks = []
	areas = []
	norm_areas = []
	for l in mrna_lens:
		profiles.append(fm.frag_factor_profile(l, frag_dist=frag_dist, no_end=no_end))
		# Calculate the peaks and areas
		peaks.append(profiles[-1].max())
		areas.append(profiles[-1].sum())
		norm_areas.append(profiles[-1].sum()/len(profiles[-1]))
	return profiles, peaks, areas, norm_areas

###############################################################################

DATA_PREFIX = '../data/'
mrna_lens = [1]+range(50,1200,40)

# Make sure we are using most up to date files
resave_frag_dists(DATA_PREFIX + 'seq_to_cluster_used.csv', 
	              DATA_PREFIX + 'stats_raw/', 
	              DATA_PREFIX + 'stats/')

# Plot the results
fig = plt.figure(figsize=(5,2.7))
# Color and alpha for the lines
c = (0.0,0.0,0.0)
a = 0.03
# Set up the first axis
ax1 = plt.subplot(1,2,1)
ax1.set_title('', fontsize=12)
ax1.set_xlabel('mRNA Length (bp)', fontsize=12)
ax1.set_ylabel('Max fraction of actual reads', fontsize=12)
ax1.tick_params(axis='x', labelsize=10)
ax1.tick_params(axis='y', labelsize=10)
plt.xlim([0,max(mrna_lens)])
plt.ylim([0,1.1])
# Set up the second axis
ax2 = plt.subplot(1,2,2)
ax2.set_title('', fontsize=12)
ax2.set_xlabel('mRNA Length (bp)', fontsize=12)
ax2.set_ylabel('Fraction of actual reads', fontsize=12)
ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)
plt.xlim([0,max(mrna_lens)])
plt.ylim([0,1.1])

# Process each variant
peaks_dist = []
areas_dist = []
norm_area_dist = []
for i in range(1,85):
	# Replicate 1
	print 'Processing variant:', i, '(rep 1)'
	reads, reads_norm = read_frag_dist('./data/stats/variant_'+str(i)+'_1.stats')
	profiles, peaks, areas, norm_areas = generate_frag_data(reads_norm, 
	                                                        mrna_lens)
	# Plot the data
	ax1.plot(mrna_lens, peaks, color=c, alpha=a)
	ax2.plot(mrna_lens, norm_areas, color=c, alpha=a)
	peaks_dist.append(peaks)
	areas_dist.append(areas)
	norm_area_dist.append(norm_areas)
	# Replicate 2
	print 'Processing variant:', i, '(rep 2)'
	reads, reads_norm = read_frag_dist('./data/stats/variant_'+str(i)+'_2.stats')
	profiles, peaks, areas, norm_areas = generate_frag_data(reads_norm, 
	                                                        mrna_lens)
	# Plot the data
	ax1.plot(mrna_lens, peaks, color=c, alpha=a)
	ax2.plot(mrna_lens, norm_areas, color=c, alpha=a)
	peaks_dist.append(peaks)
	areas_dist.append(areas)
	norm_area_dist.append(norm_areas)

# Save the plot
plt.tight_layout()
plt.savefig('testing.pdf')


