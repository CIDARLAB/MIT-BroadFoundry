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

###############################################################################
# Nif specific analysis
###############################################################################

nif_names = ['nifW', 'nifZ', 'nifM', 'nifU', 'nifV', 'nifS']
nif_lens = [10, 20, 50, 100, 200, 300]

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
plt.tight_layout()
plt.savefig('test.pdf')

# Clear the plotting cache
plt.close('all')
