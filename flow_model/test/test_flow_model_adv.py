#!/usr/bin/env python
"""
Tests of the flow model (nonhomogenious promoter and terminator rates)
"""
__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import gene_cluster_library as gcl
import gene_cluster_visualization as gcv
import flow_model as flow

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# Load the Stata nif library data
nifs = gcl.GeneClusterLibrary()
nifs.load('../../data_sets/clean_nif_stata_library.txt')

# Variant to simulate
variant = '75'

# Test the model
model = flow.generate_nonhomogeneous_site_model(nifs, variant, 1, -2, site_len=25,
	                                            rates=flow.FLOW_MODEL_RATES, 
	                                            p_rate_attrib='REU', p_rate_factor=0.5,
	                                            t_rate_attrib='Strength', t_rate_factor=0.2)

# Run the simulation
y, info = flow.run_flow_model(model, converged_site_err=0.005, sim_step_time=20.0,
	                     max_sim_time=10000.0, min_iter=20, verbose=True)

# Plot the results with architecture (including RNA-seq data)
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
fig = plt.figure(figsize=(14,4))
ax_arch = plt.subplot(gs[1])
ax_traces = plt.subplot(gs[0])

# Load the rnap densities and plot
ts = []
rnap_len = len(model[0][0])
ts.append(y[0:rnap_len])
ts.append(y[rnap_len:])
trace_len = len(ts[0])
ax_traces.fill_between(range(trace_len),ts[0],np.zeros(trace_len), color='pink', 
	                   edgecolor='red', linewidth=1.2, zorder=1)
ax_traces.fill_between(range(trace_len),-ts[1],np.zeros(trace_len), color='lightblue', 
	                   edgecolor='blue', linewidth=1.2, zorder=1)
ax_traces.plot(range(trace_len), np.zeros(trace_len), color=(0,0,0), linewidth=1.2, zorder=2)

# Scale the y-axis of the traces appropriately
max_read_depth = max(ts[0])
max_read_depth_1 = max(ts[1])
if max_read_depth_1 > max_read_depth:
	max_read_depth = max_read_depth_1
max_read_depth *= 1.05

# Update axis visibility
ax_traces.set_ylim([-max_read_depth,max_read_depth])
ax_traces.set_xlim([0,rnap_len])
ax_traces.spines["right"].set_visible(False)
ax_traces.spines["top"].set_visible(False)
ax_traces.spines["bottom"].set_visible(False)
ax_traces.tick_params(axis='both', direction='out')
ax_traces.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax_traces.set_xticks([])
ax_traces.get_yaxis().tick_left()
ax_traces.tick_params(axis='x', labelsize=8)

# Color mapping to use in library rendering
cmap = {}
# Promoters
cmap['P1'] = (1.0,0.8,0.8)
cmap['P2'] = (1.0,0.5,0.5)
cmap['P3'] = (1.0,0.0,0.0)
# Terminator
cmap['T1'] = (0.0,0.0,1.0)
# nif Genes
cmap['nifM'] = (0.5,0.77,0.56)
cmap['nifS'] = (0.37,0.78,0.78)
cmap['nifU'] = (0.76,0.47,0.76)
cmap['nifV'] = (0.78,0.78,0.49)
cmap['nifW'] = (0.94,0.74,0.4)
cmap['nifZ'] = (0.91,0.6,0.58)

# Plot the architecture below (scaling should be correct)
gcv.plot_variant_arch(ax_arch, nifs, variant, start_idx=1, end_idx=-2, linewidth=1.2, colormap=cmap)
plt.tight_layout()
fig.savefig('flow_example.pdf')
plt.show()
