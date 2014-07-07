#!/usr/bin/env python
"""
Gene Cluster Visualization
==========================

    This module contains functions to enable easier visualization of Gene Cluster
    Libraries. It can plot traces in addition to standard promoter, RBS, gene 
    lke figures. To simplify export of results, all plotting is performed by 
    matplotlib which ensures portability.
"""
#    Gene Cluster Visualization
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
nifs.load('../data/nif_stata_library.txt')

# Example plot testing promoter/CDS/terminator renderings
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(1,1,1)
ax.plot([0,200], [0,0], color=(0,0,0), linewidth=2, zorder=1)
gcv.draw_promoter(ax, 5, 25)
gcv.draw_cds(ax, 35, 75)
gcv.draw_terminator(ax, 76, 90)
gcv.draw_promoter(ax, 195, 155)
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

# Example plot of the whole library
fig = plt.figure(figsize=(14,84))
gcv.plot_library_arch(fig, nifs, linewidth=1.2)
plt.tight_layout()
fig.savefig('./visualizations/LibraryRendering.pdf')

# Example plot of architecture with Tx traces
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
fig = plt.figure(figsize=(14,4))
ax_arch = plt.subplot(gs[1])
ax_traces = plt.subplot(gs[0],sharex=ax_arch)
# Load the traces for predicted and measured
phys_reads = gca.load_strand_data('../data/phys_depths3.csv')
gcv.plot_traces_with_arch(ax_arch, ax_traces, nifs, '25', phys_reads, start_idx=1, end_idx=-1, linewidth=1.2)
plt.tight_layout()
fig.savefig('./visualizations/ArchAndTrace.pdf')

# Clear the plotting cache
plt.close('all')
