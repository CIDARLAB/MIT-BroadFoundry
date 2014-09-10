#!/usr/bin/env python
"""
Analyse the fragmentation model
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import fragmentation_model as fm
import matplotlib.pyplot as plt
import numpy as np

# Statistical distribution of fragments
frag_mean = 280
frag_sd = 70

GENERAL_ANALYSIS = True

if GENERAL_ANALYSIS == True:
	# Generate range of profiles for differing length mRNAs
	mrna_lens = [1] + range(50,1050,50) + [1050, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 2000]
	profiles = []
	peaks = []
	areas = []
	for l in mrna_lens:
		print 'Generating profile for mRNA of length:', l
		profiles.append(fm.frag_factor_profile(l, frag_mean, frag_sd, bp_cutoff=10*frag_sd, no_end=False))
		# Calculate the peaks and areas
		peaks.append(profiles[-1].max())
		areas.append(profiles[-1].sum())

	# Plot the results
	fig = plt.figure(figsize=(8,4))

	ax1 = plt.subplot(1,2,1)
	ax1.plot(mrna_lens, peaks)
	plt.title('Read Depth Peak', fontsize=12)
	plt.xlabel('mRNA Length (bp)', fontsize=12)
	plt.ylabel('Read Depth Peak (a.u.)', fontsize=12)
	plt.xlim([0,max(mrna_lens)])
	#ax1.set_yscale('log')

	ax2 = plt.subplot(1,2,2)
	ax2.plot(mrna_lens, areas)
	plt.title('Read Depth Area', fontsize=12)
	plt.xlabel('mRNA Length (bp)', fontsize=12)
	plt.ylabel('Read Depth Area (a.u.)', fontsize=12)
	plt.xlim([0,max(mrna_lens)])
	#ax2.set_yscale('log')

	plt.tight_layout()
	plt.savefig('mrna_length_rnaseq_signal.pdf')

###############################################################################
# Nif specific analysis
###############################################################################

nif_names = ['nifW', 'nifZ', 'nifM', 'nifU', 'nifV', 'nifS']
nif_lens = [258, 447, 801, 825, 1143, 1203]

profiles = []
peaks = []
areas = []
for l in nif_lens:
	print 'Generating nif profile for mRNA of length:', l
	profiles.append(fm.frag_factor_profile(l, frag_mean, frag_sd, bp_cutoff=10*frag_sd, no_end=False))
	# Calculate the peaks and areas
	peaks.append(profiles[-1].max())
	areas.append(profiles[-1].sum())

# Plot the results
fig = plt.figure(figsize=(6,4))
ax = plt.subplot(1,1,1)
for p in profiles:
	ax.plot(range(len(p)), p)
plt.title('Read Depth Peak', fontsize=12)
plt.xlabel('mRNA Length (bp)', fontsize=12)
plt.ylabel('Read Depth Peak (a.u.)', fontsize=12)
plt.xlim([0,1203])
plt.tight_layout()
plt.savefig('nif_length_analysis.pdf')

# Clear the plotting cache
plt.close('all')
