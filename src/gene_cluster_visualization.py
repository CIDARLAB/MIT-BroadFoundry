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

import sys
if sys.version_info[:2] < (2, 6):
    m = "Python version 2.6 or later is required for Gene Cluster Tools (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

import csv
import numpy as np
import matplotlib.pyplot as plt
import gene_cluster_library as gcl

from matplotlib.patches import Polygon
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def draw_promoter (ax, start_bp, end_bp, color=(1.0,0.0,0.0), extent=50, offset=10, linewidth=2, 
	               p_width=1, arrow_len=6, arrow_width=1, label=None):
	# Draw the promoter symbol
	dir_fac = 1.0
	if start_bp > end_bp:
		dir_fac = -1.0
	l1 = Line2D([start_bp,start_bp],[0,dir_fac*offset], linewidth=linewidth, color=color, zorder=10)
	l2 = Line2D([start_bp,start_bp+dir_fac*extent-dir_fac*(arrow_len*0.5)],[dir_fac*offset,dir_fac*offset], linewidth=linewidth, color=color, zorder=10)
	ax.add_line(l1)
	ax.add_line(l2)
	p1 = Polygon([(start_bp+dir_fac*extent-dir_fac*arrow_len, dir_fac*offset+(arrow_width)), 
		          (start_bp+dir_fac*extent, dir_fac*offset),
		          (start_bp+dir_fac*extent-dir_fac*arrow_len, dir_fac*offset-(arrow_width))],
		          facecolor=color, edgecolor=color, linewidth=linewidth)
	ax.add_patch(p1)
 	# Shade the promoter area (normally smaller than symbol extent)
 	p2 = Polygon([(start_bp, -p_width), 
 		          (start_bp, p_width),
 		          (end_bp, p_width),
 		          (end_bp, -p_width)], facecolor=color, edgecolor=color, linewidth=linewidth, zorder=10)
	ax.add_patch(p2)

def draw_terminator (ax, start_bp, end_bp, color=(0.0,0.0,1.0), extent=5, offset=6, linewidth=2, t_width=1, label=None):
	# Draw the terminator symbol
	dir_fac = 1.0
	if start_bp > end_bp:
		dir_fac = -1.0
	l1 = Line2D([start_bp,start_bp],[0,dir_fac*offset], linewidth=linewidth, color=color, zorder=9)
	l2 = Line2D([start_bp-extent,start_bp+extent],[dir_fac*offset,dir_fac*offset], linewidth=linewidth, color=color, zorder=9)
	ax.add_line(l1)
	ax.add_line(l2)
	# Shade the terminator area (normally smaller than symbol extent)
 	p2 = Polygon([(start_bp, -t_width), 
 		          (start_bp, t_width),
 		          (end_bp, t_width),
 		          (end_bp, -t_width)], facecolor=color, edgecolor=color, linewidth=linewidth, zorder=9)
	ax.add_patch(p2)

def draw_cds (ax, start_bp, end_bp, color=(0.0,1.0,0.0), linewidth=2, arrow_extent=16, body_offset=3, arrow_offset=3, 
	          hatch='', label=None):
	# Draw the CDS symbol
	dir_fac = 1.0
	if start_bp > end_bp:
		dir_fac = -1.0
	p1 = Polygon([(start_bp, body_offset), 
		          (start_bp, -body_offset),
		          (end_bp-dir_fac*arrow_extent, -body_offset),
		          (end_bp-dir_fac*arrow_extent, -body_offset-arrow_offset),
		          (end_bp, 0),
		          (end_bp-dir_fac*arrow_extent, body_offset+arrow_offset),
		          (end_bp-dir_fac*arrow_extent, body_offset)],
		          facecolor=color, edgecolor=(0.0,0.0,0.0), linewidth=linewidth, hatch=hatch, zorder=8)
	ax.add_patch(p1)

def plot_library_arch (fig, gcl, linewidth=1.0):
	for v_idx in range(len(gcl.variants.keys())):
		v = gcl.variants.keys()[v_idx]
		ax = fig.add_subplot(len(gcl.variants.keys()),1,v_idx)
		plot_variant_arch(ax, gcl, v, start_idx=1, end_idx=-1, linewidth=linewidth)

def plot_variant_arch (ax, gcl, variant, start_idx=0, end_idx=None, linewidth=1.0):
	# Calculate start and end bp
	part_list = gcl.variants[variant]['part_list']
	if end_idx == None:
		end_idx = len(part_list)-1
	if end_idx < 0:
		end_idx = len(part_list)-1+end_idx
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
		if part_type == 'Promoter':
			draw_promoter(ax, part_start_bp, part_end_bp, linewidth=linewidth, extent=60, arrow_len=30, arrow_width=1)
		if part_type == 'CDS':
			draw_cds(ax, part_start_bp, part_end_bp, linewidth=linewidth, arrow_extent=80)
		if part_type == 'Terminator':
			draw_terminator(ax, part_start_bp, part_end_bp, linewidth=linewidth, extent=20)
	# Resize the axis and set bounds to ensure fully visible
	start_bp = part_list[start_idx]['seq_idx']
	if part_list[0]['dir'] == 'R':
		start_bp -= part_list[start_idx]['seq_len']
	end_bp = part_list[end_idx]['seq_idx'] + part_list[end_idx]['seq_len']
	if part_list[0]['dir'] == 'R':
		end_bp -= part_list[end_idx]['seq_len']
	end_padding = 20
	ax.plot([start_bp-end_padding,end_bp+end_padding], [0,0], color=(0,0,0), linewidth=linewidth, zorder=1)
	ax.set_xlim([start_bp-end_padding,end_bp+end_padding])
	ax.set_ylim([-15,15])
	ax.set_axis_off()

def plot_traces_with_arch (ax_arch, ax_traces, gcl, variant, traces, start_idx=None, end_idx=None, linewidth=1.2):
	# Plot the architecture of the cluster
	plot_variant_arch(ax_arch, gcl, variant, start_idx=start_idx, end_idx=end_idx, linewidth=linewidth)
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

	ax_traces.set_ylim([-max_read_depth,max_read_depth])

	ax_traces.spines["right"].set_visible(False)
	ax_traces.spines["top"].set_visible(False)
	ax_traces.spines["bottom"].set_visible(False)
	ax_traces.tick_params(axis='both', direction='out')
	ax_traces.get_xaxis().tick_bottom()   # remove unneeded ticks 
	ax_traces.set_xticks([])
	ax_traces.get_yaxis().tick_left()
	ax_traces.tick_params(axis='x', labelsize=8)

def load_strand_data (filename):
	data = {}
	reader = csv.reader(open(filename, 'rb'), delimiter=',')
	# Use the header to get the names of the keys for the prediction data
	header = reader.next()
	for el in header[1:]:
		el_name = el[6:-1]
		# Only consider the first of the replicates
		if el_name.count('.') == 1:
			el_name = el_name[0]
		if el_name not in data.keys():
			data[el_name] = [[],[]]
	for row in reader:
		for el_idx in range(1, len(row)):
			el_name = header[el_idx][6:-1]
			# Only consider the first of the replicates
			if el_name.count('.') == 1:
				el_name = el_name[0]
			strand = header[el_idx][-1]
			if strand == '+':
				data[el_name][0].append(float(row[el_idx]))
			else:
				data[el_name][1].append(float(row[el_idx]))
	return data

# Load the Stata nif library data
nifs = gcl.GeneClusterLibrary()
nifs.load('./data/nif_stata_library.txt')

# Example plot testing promoter/CDS/terminator renderings
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(1,1,1)
ax.plot([0,200], [0,0], color=(0,0,0), linewidth=2, zorder=1)
draw_promoter(ax, 5, 25)
draw_cds(ax, 35, 75)
draw_terminator(ax, 76, 90)
draw_promoter(ax, 195, 155)
draw_cds(ax, 150, 120, hatch='////')
draw_terminator(ax, 119, 100)
ax.set_xlim([0,200])
ax.set_ylim([-30,30])
ax.set_axis_off()

# Example plot testing full variant rendering
fig = plt.figure(figsize=(14,3))
ax = fig.add_subplot(3,1,1)
plot_variant_arch(ax, nifs, '25', start_idx=1, end_idx=-1, linewidth=1.2)
ax = fig.add_subplot(3,1,2)
plot_variant_arch(ax, nifs, '10', start_idx=1, end_idx=-1, linewidth=1.2)
ax = fig.add_subplot(3,1,3)
plot_variant_arch(ax, nifs, '1', start_idx=1, end_idx=-1, linewidth=1.2)
plt.tight_layout()

# Example plot of the whole library
fig = plt.figure(figsize=(14,84))
plot_library_arch(fig, nifs, linewidth=1.2)
plt.tight_layout()

# Example plot of architecture with Tx traces
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
fig = plt.figure(figsize=(14,4))
ax_arch = plt.subplot(gs[1])
ax_traces = plt.subplot(gs[0],sharex=ax_arch)
# Load the traces for predicted and measured
phys_reads = load_strand_data('./data/phys_depths3.csv')
plot_traces_with_arch(ax_arch, ax_traces, nifs, '10', phys_reads, start_idx=1, end_idx=-1, linewidth=1.2)
plt.tight_layout()

# Show the examples
plt.show()

