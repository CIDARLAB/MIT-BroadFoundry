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
frag_mean = 280
frag_sd = 70

# Generate range of profiles for differing length mRNAs
mrna_lens = [1] + range(50,1200,40) # 1200
profiles = []
peaks = []
areas = []
norm_areas = []
for l in mrna_lens:
	print 'Generating profile for mRNA of length:', l
	profiles.append(fm.frag_factor_profile(l, frag_mean, frag_sd, bp_cutoff=10*frag_sd, no_end=False))
	# Calculate the peaks and areas
	peaks.append(profiles[-1].max())
	areas.append(profiles[-1].sum())
	norm_areas.append(profiles[-1].sum()/len(profiles[-1]))

# Plot the results
fig = plt.figure(figsize=(4,2.5))

ax1 = fig.add_axes([0.145,0.2,0.7,0.75], frameon=True)
#ax1 = plt.subplot(1,1,1)
my_cmap = cm.get_cmap('RdYlBu')
N = len(profiles)
c = 1
for i in xrange(N):
	c = my_cmap(float(i)/(N-1))
	ax1.plot(range(len(profiles[i])), profiles[i], color=c)
plt.title('', fontsize=12)
plt.xlabel('mRNA Position (bp)', fontsize=12)
plt.ylabel('Fraction of actual reads', fontsize=12)
ax1.tick_params(axis='x', labelsize=10)
ax1.tick_params(axis='y', labelsize=10)
plt.xlim([0,max(mrna_lens)])
plt.ylim([0,1.1])
# Plot colorbar.
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.normalize(vmin=0, vmax=1))
# fake up the array of the scalar mappable. Urgh...
sm._A = []
axcolor = fig.add_axes([0.865,0.2,0.03,0.75])
cbar1 = pylab.colorbar(sm, cax=axcolor)
cbar1.ax.tick_params(labelsize=10)
tick_locs = [0.0, 1.0]
tick_labels = ['1','1200']
cbar1.locator = matplotlib.ticker.FixedLocator(tick_locs)
cbar1.formatter = matplotlib.ticker.FixedFormatter(tick_labels)
cbar1.update_ticks()
cbar1.ax.set_ylabel('mRNA length (bp)', fontsize=10, rotation=270, labelpad=-21)

plt.savefig('mrna_length_rnaseq_signal_2.pdf')

# Plot the results
fig = plt.figure(figsize=(5,2.7))

ax2 = plt.subplot(1,2,1)
ax2.plot(mrna_lens, peaks, color=c)
plt.title('', fontsize=12)
plt.xlabel('mRNA Length (bp)', fontsize=12)
plt.ylabel('Max fraction of actual reads', fontsize=12)
ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)
plt.xlim([0,max(mrna_lens)])
plt.ylim([0,1.1])

ax3 = plt.subplot(1,2,2)
ax3.plot(mrna_lens, norm_areas, color=c)
plt.title('', fontsize=12)
plt.xlabel('mRNA Length (bp)', fontsize=12)
plt.ylabel('Fraction of actual reads', fontsize=12)
ax3.tick_params(axis='x', labelsize=10)
ax3.tick_params(axis='y', labelsize=10)
plt.xlim([0,max(mrna_lens)])
plt.ylim([0,1.1])

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
	profiles.append(fm.frag_factor_profile(l, frag_mean, frag_sd, bp_cutoff=None, no_end=False))
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
plt.ylim([0,1.1])
plt.tight_layout()
plt.savefig('nif_length_analysis.pdf')

# Clear the plotting cache
plt.close('all')
