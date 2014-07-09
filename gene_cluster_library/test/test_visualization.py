#!/usr/bin/env python
"""
Test the visualization functionality of the library.
"""
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gene_cluster_library as gcl
import gene_cluster_analysis as gca
import gene_cluster_visualization as gcv

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

# Load the Stata nif library data
nifs = gcl.GeneClusterLibrary()
nifs.load('./data/nif_stata_library.txt')

# Example plot testing promoter/CDS/terminator renderings
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(1,1,1)
ax.plot([0,200], [0,0], color=(0,0,0), linewidth=2, zorder=1)
gcv.draw_promoter(ax, 5, 25)
gcv.draw_rbs(ax, 29, 34)
gcv.draw_cds(ax, 35, 75)
gcv.draw_terminator(ax, 76, 90)
gcv.draw_promoter(ax, 195, 155)
gcv.draw_rbs(ax, 154, 151)
gcv.draw_cds(ax, 150, 120, hatch='////')
gcv.draw_terminator(ax, 119, 100)
ax.set_xlim([0,200])
ax.set_ylim([-30,30])
ax.set_axis_off()
fig.savefig('./visualizations/PartDrawing.pdf')

# Example plot testing full variant rendering
fig = plt.figure(figsize=(14,3))
ax = fig.add_subplot(3,1,1)
gcv.plot_variant_arch(ax, nifs, '25', start_idx=1, end_idx=-1, linewidth=1.2)
ax = fig.add_subplot(3,1,2)
gcv.plot_variant_arch(ax, nifs, '10', start_idx=1, end_idx=-1, linewidth=1.2)
ax = fig.add_subplot(3,1,3)
gcv.plot_variant_arch(ax, nifs, '1', start_idx=1, end_idx=-1, linewidth=1.2)
plt.tight_layout()
fig.savefig('./visualizations/VariantRendering.pdf')

# Color mapping to use in library rendering
cmap = {}
# Promoters
cmap['P1'] = (1.0,0.0,0.0)
cmap['P2'] = (1.0,0.5,0.5)
cmap['P3'] = (1.0,0.8,0.8)
# Terminator
cmap['T1'] = (0.0,0.0,1.0)
# nif Genes
cmap['nifM'] = (0.5,0.77,0.56)
cmap['nifS'] = (0.37,0.78,0.78)
cmap['nifU'] = (0.76,0.47,0.76)
cmap['nifV'] = (0.78,0.78,0.49)
cmap['nifW'] = (0.94,0.74,0.4)
cmap['nifZ'] = (0.91,0.6,0.58)

rbsStrong = (0.16,0.68,0.15)
rbsWeak = (0.47,0.83,0.46)
cmap['Rm1'] = rbsStrong
cmap['Rm2'] = rbsWeak
cmap['Rs1'] = rbsStrong
cmap['Rs2'] = rbsWeak
cmap['Ru1'] = rbsStrong
cmap['Ru2'] = rbsWeak
cmap['Rv1'] = rbsStrong
cmap['Rv2'] = rbsWeak
cmap['Rw1'] = rbsStrong
cmap['Rw2'] = rbsWeak
cmap['Rz1'] = rbsStrong
cmap['Rz2'] = rbsWeak

# Hatch mapping to use in library rendering
hmap = {}
hmap['nifM'] = ''
hmap['nifS'] = '/'
hmap['nifU'] = '//'
hmap['nifV'] = '///'
hmap['nifW'] = '////'
hmap['nifZ'] = '/////'

# Example plot of the whole library
fig = plt.figure(figsize=(14,84))
gcv.plot_library_arch(fig, nifs, linewidth=1.2, colormap=cmap, hatchmap=hmap)
plt.tight_layout()
fig.savefig('./visualizations/LibraryRendering.pdf')

# Example plot of architecture with Tx traces
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
fig = plt.figure(figsize=(14,4))
ax_arch = plt.subplot(gs[1])
ax_traces = plt.subplot(gs[0],sharex=ax_arch)
# Load the traces for predicted and measured
phys_reads = gca.load_strand_data('./data/phys_depths3.csv')
gcv.plot_traces_with_arch(ax_arch, ax_traces, nifs, '25', phys_reads, start_idx=1, end_idx=-1, linewidth=1.2, colormap=cmap, hatchmap=hmap)
plt.tight_layout()
fig.savefig('./visualizations/ArchAndTrace.pdf')

# Clear the plotting cache
plt.close('all')
