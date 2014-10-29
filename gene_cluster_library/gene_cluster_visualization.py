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

def plot_library_arch (fig, gcl, start_idx=0, end_idx=None, linewidth=1.0, scaleadjust=1.0):
	for v_idx in range(len(gcl.variants.keys())):
		v = gcl.variants.keys()[v_idx]
		ax = fig.add_subplot(len(gcl.variants.keys()),1,v_idx)
		plot_variant_arch(ax, gcl, v, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth, scaleadjust=scaleadjust)

def plot_variant_arch (ax, gcl, variant, start_idx=0, end_idx=None, linewidth=1.0, end_padding=20, scaleadjust=1.0):
	# Calculate start and end bp
	part_list = gcl.variants[variant]['part_list']
	if end_idx == None:
		end_idx = len(part_list)-1
	if end_idx < 0:
		end_idx = len(part_list)+end_idx
	# Build the object for rendering
	design = []
	for el_idx in range(start_idx, end_idx+1):
		part_name = part_list[el_idx]['part_name']
		part_type =  gcl.parts[part_name]['type']
		part_start_bp = part_list[el_idx]['seq_idx']
		part_end_bp = part_start_bp
		if part_list[el_idx]['dir'] == 'F':	
			part_end_bp += part_list[el_idx]['seq_len']
		else:
			part_start_bp += part_list[el_idx]['seq_len']
		opts = {}
		for k in part_list[el_idx].keys():
			if k not in ['part_name', 'type', 'seq_idx', 'dir']:
				opts[k] = part_list[el_idx][k]
		design.append({'type':part_type, 'name':part_name, 'fwd':part_list[el_idx]['dir'], 'start':part_start_bp, 'end':part_end_bp, 'opts':opts})
	# Find the min and max bp to scale things properly
	min_bp = design[0]['start']
	max_bp = design[0]['start']
	for el in design:
		if el['start'] < min_bp:
			min_bp = el['start']
		if el['end'] < min_bp:
			min_bp = el['end']
		if el['start'] > max_bp:
			max_bp = el['start']
		if el['end'] > max_bp:
			max_bp = el['end']
	# Draw the elements
	design_len = float(max_bp-min_bp)
	# Create the DNA renderer
	dr = dpl.DNARenderer(scale=(design_len/2000.0)*scaleadjust, linewidth=linewidth)
	# Set the renders to draw elements (we don't draw regulation)
	part_renderers = dr.trace_part_renderers()
	# Redender the DNA to axis
	start, end = dr.renderDNA(ax, design, part_renderers)
	# Resize the axis and set bounds to ensure fully visible
	start_bp = part_list[start_idx]['seq_idx']
	if part_list[0]['dir'] == 'R':
		start_bp -= part_list[start_idx]['seq_len']
	end_bp = part_list[end_idx]['seq_idx'] + part_list[end_idx]['seq_len']
	if part_list[0]['dir'] == 'R':
		end_bp -= part_list[end_idx]['seq_len']
	# Plot the backbone
	bb_min_bp = design[0]['start']
	bb_max_bp = design[0]['start']
	for el in part_list:
		temp_start_bp = el['seq_idx']
		temp_end_bp = temp_start_bp
		if el['dir'] == 'F':	
			temp_end_bp += el['seq_len']
		else:
			temp_start_bp += el['seq_len']
		if temp_start_bp < bb_min_bp:
			bb_min_bp = temp_start_bp
		if temp_end_bp < bb_min_bp:
			bb_min_bp = temp_end_bp
		if temp_start_bp > bb_max_bp:
			bb_max_bp = temp_start_bp
		if temp_end_bp > bb_max_bp:
			bb_max_bp = temp_end_bp
	ax.plot([bb_min_bp,bb_max_bp], [0,0], color=(0,0,0), linewidth=linewidth, zorder=1)
	ax.set_xlim([start_bp-end_padding,end_bp+end_padding])
	ax.set_ylim([-8,8])
	ax.set_axis_off()

def plot_traces_with_arch (ax_arch, ax_traces, gcl, variant, traces, start_idx=None, end_idx=None, linewidth=1.2, scaleadjust=1.0):
	# Plot the architecture of the cluster
	plot_variant_arch(ax_arch, gcl, variant, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth, scaleadjust=scaleadjust)
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
	ax_traces.get_xaxis().tick_bottom()
	ax_traces.set_xticks([])
	ax_traces.get_yaxis().tick_left()
	ax_traces.tick_params(axis='x', labelsize=8)

def plot_dual_traces_with_arch (ax_arch, ax_traces_1, ax_traces_2, gcl, variant, traces_1, traces_2, start_idx=None, end_idx=None, linewidth=1.2, scaleadjust=1.0):
	# Plot the architecture of the cluster
	plot_variant_arch(ax_arch, gcl, variant, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth, scaleadjust=scaleadjust)
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
	ax_traces_2.get_xaxis().tick_bottom()
	ax_traces_2.set_xticks([])
	ax_traces_2.get_yaxis().tick_left()
	ax_traces_2.tick_params(axis='x', labelsize=8)

