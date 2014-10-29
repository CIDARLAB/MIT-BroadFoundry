#!/usr/bin/env python
"""
Gene Cluster Visualization
==========================

    This module contains functions to enable easier visualization of Gene Cluster
    Libraries. It can plot traces in addition to standard promoter, gene, terminator 
    like figures. To simplify export of results, all plotting is performed by 
    matplotlib which ensures portability and vector based output for publication
    quality figures.
"""
#    Gene Cluster Visualization
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import sys
if sys.version_info[:2] < (2, 6):
    m = "Python version 2.6 or later is required for Gene Cluster Tools (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

import numpy as np
import gene_cluster_library as gcl
import dnaplotlib as dpl

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def plot_library_arch (fig, gcl, start_idx=0, end_idx=None, linewidth=1.0, colormap=None, hatchmap=None):
	for v_idx in range(len(gcl.variants.keys())):
		v = gcl.variants.keys()[v_idx]
		ax = fig.add_subplot(len(gcl.variants.keys()),1,v_idx)
		plot_variant_arch(ax, gcl, v, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth, colormap=colormap, hatchmap=hatchmap)

def plot_variant_arch (ax, gcl, variant, start_idx=0, end_idx=None, linewidth=1.0, colormap=None, hatchmap=None):
	# Calculate start and end bp
	part_list = gcl.variants[variant]['part_list']
	if end_idx == None:
		end_idx = len(part_list)-1
	if end_idx < 0:
		end_idx = len(part_list)+end_idx
	# Draw the parts
	for el_idx in range(start_idx, end_idx+1):
		part_name = part_list[el_idx]['part_name']
		part_type =  gcl.parts[part_name]['type']
		part_start_bp = part_list[el_idx]['seq_idx']
		part_end_bp = part_start_bp
		if part_list[el_idx]['dir'] == 'F':	
			part_end_bp += part_list[el_idx]['seq_len']
		else:
			part_start_bp += part_list[el_idx]['seq_len']
		part_color = None
		if colormap != None and part_name in colormap.keys():
			part_color = colormap[part_name]
		part_hatch = None
		if hatchmap != None and part_name in hatchmap.keys():
			part_hatch = hatchmap[part_name]
		
		# Build the object for rendering

		




		if part_type == 'Promoter':
			draw_promoter(ax, part_start_bp, part_end_bp, color=part_color, linewidth=linewidth, extent=60, arrow_len=30, arrow_width=1)
		if part_type == 'CDS':
			draw_cds(ax, part_start_bp, part_end_bp, color=part_color, hatch=part_hatch, linewidth=linewidth, arrow_extent=80)
		if part_type == 'Terminator':
			draw_terminator(ax, part_start_bp, part_end_bp, color=part_color, linewidth=linewidth, extent=20)
		if part_type == 'RBS':
			draw_rbs(ax, part_start_bp, part_end_bp, color=part_color, linewidth=linewidth, extent=30)
	



	# Draw the elements
	scale = 

	# Create the DNA renderer
	dr = dpl.DNARenderer(scale=2000.0/1500.0) # Recomend scaling bp by 1500.0 to get similar glyph sizes

	# Set the renders to draw elements
	part_renderers = dr.trace_part_renderers()


	# Resize the axis and set bounds to ensure fully visible
	start_bp = part_list[start_idx]['seq_idx']
	if part_list[0]['dir'] == 'R':
		start_bp -= part_list[start_idx]['seq_len']
	end_bp = part_list[end_idx]['seq_idx'] + part_list[end_idx]['seq_len']
	if part_list[0]['dir'] == 'R':
		end_bp -= part_list[end_idx]['seq_len']
	end_padding = 20
	#ax.plot([start_bp-end_padding,end_bp+end_padding], [0,0], color=(0,0,0), linewidth=linewidth, zorder=1)
	ax.set_xlim([start_bp-end_padding,end_bp+end_padding])
	ax.set_ylim([-15,15])
	ax.set_axis_off()

def plot_traces_with_arch (ax_arch, ax_traces, gcl, variant, traces, start_idx=None, end_idx=None, linewidth=1.2, colormap=None, hatchmap=None):
	# Plot the architecture of the cluster
	plot_variant_arch(ax_arch, gcl, variant, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth, colormap=colormap, hatchmap=hatchmap)
	# Plot the traces (shared axis ensures only valid region plotted)
	ts = np.array(traces[variant])
	trace_len = len(ts[0])
	ax_traces.fill_between(range(trace_len),ts[0],np.zeros(trace_len), color='pink', edgecolor='red', linewidth=linewidth, zorder=1)
	ax_traces.fill_between(range(trace_len),-ts[1],np.zeros(trace_len), color='lightblue', edgecolor='blue', linewidth=linewidth, zorder=1)
	ax_traces.plot(range(trace_len), np.zeros(trace_len), color=(0,0,0), linewidth=linewidth, zorder=2)
	# Scale the y-axis of the traces appropriately
	max_read_depth = max(ts[0])
	max_read_depth_1 = max(ts[1])
	if max_read_depth_1 > max_read_depth:
		max_read_depth = max_read_depth_1
	max_read_depth *= 1.05
	# Update axis visibility
	ax_traces.set_ylim([-max_read_depth,max_read_depth])
	#ax_traces.set_xlim([0,trace_len])
	ax_traces.spines["right"].set_visible(False)
	ax_traces.spines["top"].set_visible(False)
	ax_traces.spines["bottom"].set_visible(False)
	ax_traces.tick_params(axis='both', direction='out')
	ax_traces.get_xaxis().tick_bottom()   # remove unneeded ticks 
	ax_traces.set_xticks([])
	ax_traces.get_yaxis().tick_left()
	ax_traces.tick_params(axis='x', labelsize=8)

def plot_dual_traces_with_arch (ax_arch, ax_traces_1, ax_traces_2, gcl, variant, traces_1, traces_2, start_idx=None, end_idx=None, linewidth=1.2, colormap=None, hatchmap=None):
	# Plot the architecture of the cluster
	plot_variant_arch(ax_arch, gcl, variant, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth, colormap=colormap, hatchmap=hatchmap)
	# Plot the traces (shared axis ensures only valid region plotted)
	ts = np.array(traces_1[variant])
	trace_len = len(ts[0])
	ax_traces_1.fill_between(range(trace_len),ts[0],np.zeros(trace_len), color='pink', edgecolor='red', linewidth=linewidth, zorder=1)
	ax_traces_1.fill_between(range(trace_len),-ts[1],np.zeros(trace_len), color='lightblue', edgecolor='blue', linewidth=linewidth, zorder=1)
	ax_traces_1.plot(range(trace_len), np.zeros(trace_len), color=(0,0,0), linewidth=linewidth, zorder=2)
	# Scale the y-axis of the traces appropriately
	max_read_depth = max(ts[0])
	max_read_depth_1 = max(ts[1])
	if max_read_depth_1 > max_read_depth:
		max_read_depth = max_read_depth_1
	max_read_depth *= 1.05
	# Update axis visibility
	ax_traces_1.set_ylim([-max_read_depth,max_read_depth])
	#ax_traces.set_xlim([0,trace_len])
	ax_traces_1.spines["right"].set_visible(False)
	ax_traces_1.spines["top"].set_visible(False)
	ax_traces_1.spines["bottom"].set_visible(False)
	ax_traces_1.tick_params(axis='both', direction='out')
	ax_traces_1.get_xaxis().tick_bottom()   # remove unneeded ticks 
	ax_traces_1.set_xticks([])
	ax_traces_1.get_yaxis().tick_left()
	ax_traces_1.tick_params(axis='x', labelsize=8)
	# Plot the second traces (shared axis ensures only valid region plotted)
	ts = np.array(traces_2[variant])
	trace_len = len(ts[0])
	ax_traces_2.fill_between(range(trace_len),ts[0],np.zeros(trace_len), color='pink', edgecolor='red', linewidth=linewidth, zorder=1)
	ax_traces_2.fill_between(range(trace_len),-ts[1],np.zeros(trace_len), color='lightblue', edgecolor='blue', linewidth=linewidth, zorder=1)
	ax_traces_2.plot(range(trace_len), np.zeros(trace_len), color=(0,0,0), linewidth=linewidth, zorder=2)
	# Scale the y-axis of the traces appropriately
	max_read_depth = max(ts[0])
	max_read_depth_1 = max(ts[1])
	if max_read_depth_1 > max_read_depth:
		max_read_depth = max_read_depth_1
	max_read_depth *= 1.05
	# Update axis visibility
	ax_traces_2.set_ylim([-max_read_depth,max_read_depth])
	#ax_traces.set_xlim([0,trace_len])
	ax_traces_2.spines["right"].set_visible(False)
	ax_traces_2.spines["top"].set_visible(False)
	ax_traces_2.spines["bottom"].set_visible(False)
	ax_traces_2.tick_params(axis='both', direction='out')
	ax_traces_2.get_xaxis().tick_bottom()   # remove unneeded ticks 
	ax_traces_2.set_xticks([])
	ax_traces_2.get_yaxis().tick_left()
	ax_traces_2.tick_params(axis='x', labelsize=8)

