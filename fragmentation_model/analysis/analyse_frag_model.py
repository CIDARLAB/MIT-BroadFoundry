#!/usr/bin/env python
"""
Analyse the fragmentation model
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import fragmentation_model as fm
import matplotlib.pyplot as plt
import pylab
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib.ticker as ticker

# Statistical distribution of fragments
frag_means = [50, 250, 500]
frag_sds   = [20, 80,  200]

# From experiments
# frag_mean = 280
# frag_sd = 70

def plot_pred_seq_coverage (ax, frag_mean, frag_sd):
	# Generate range of profiles for differing length mRNAs
	mrna_lens = range(50,1200,50)
	profiles = []
	peaks = []
	areas = []
	norm_areas = []
	for l in mrna_lens:
		print 'Generating profile for mRNA of length:', l
		profiles.append(fm.frag_factor_profile(l, frag_mean=frag_mean, frag_sd=frag_sd, no_end=False))
		# Calculate the peaks and areas
		peaks.append(profiles[-1].max())
		areas.append(profiles[-1].sum())
		norm_areas.append(profiles[-1].sum()/len(profiles[-1]))
	#ax1 = plt.subplot(1,1,1)
	my_cmap = cm.get_cmap('RdYlBu')
	N = len(profiles)
	c = 1
	plt.axhline(y=1.0, color=(0.8,0.8,0.8), linestyle='-', linewidth=1.5)
	for i in xrange(N):
		c = my_cmap(float(i)/(N-1))
		ax.plot(range(len(profiles[i])), profiles[i], color=c, linewidth=1.5)
	#plt.title('', fontsize=12)
	#plt.xlabel('mRNA Position (bp)', fontsize=12)
	#plt.ylabel('Fraction of actual reads', fontsize=12)
	ax.tick_params(axis='x', labelsize=10)
	ax.tick_params(axis='y', labelsize=10)
	plt.xlim([0,max(mrna_lens)])
	plt.ylim([0,1.19])

# Plot colorbar.
#sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.normalize(vmin=0, vmax=1))
# fake up the array of the scalar mappable. Urgh...
#sm._A = []
#axcolor = fig.add_axes([0.865,0.2,0.03,0.75])
#cbar1 = pylab.colorbar(sm, cax=axcolor)
#cbar1.ax.tick_params(labelsize=10)
#tick_locs = [0.0, 1.0]
#tick_labels = ['1','1200']
#cbar1.locator = matplotlib.ticker.FixedLocator(tick_locs)
#cbar1.formatter = matplotlib.ticker.FixedFormatter(tick_labels)
#cbar1.update_ticks()
#cbar1.ax.set_ylabel('mRNA length (bp)', fontsize=10, rotation=270, labelpad=-21)

# Plot the results
fig = plt.figure(figsize=(8.5,7))
ax_i = 1
x = 1
y = 1
for frag_mean in frag_means:
	y=1
	for frag_sd in frag_sds:
		print 'PLOTTING SUB-FIGURE:', ax_i
		ax = plt.subplot(len(frag_sds),len(frag_means),ax_i)
		plot_pred_seq_coverage(ax, frag_mean, frag_sd)
		if y != 1:
			ax.set_yticklabels([])
		if x != 3:
			ax.set_xticklabels([])
		plt.subplots_adjust(hspace=.0)
		plt.subplots_adjust(wspace=.0)
		ax_i += 1
		y += 1
	x += 1


plt.savefig('mrna_length_rnaseq_signal_2.pdf')



