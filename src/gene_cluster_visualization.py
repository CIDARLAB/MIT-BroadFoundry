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

import numpy as np
import matplotlib.pyplot as plt
import gene_cluster_library as gcl

from matplotlib.patches import Polygon
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def draw_promoter (ax, start_bp, end_bp, color=(1.0,0.0,0.0), extent=50, offset=10, linewidth=2, 
	               p_width=1, arrow_len=6, arrow_width=1, label=None):
	# Draw the promoter symbol
	dir_fac = 1.0
	if start_bp > end_bp:
		dir_fac = -1.0
	l1 = Line2D([start_bp,start_bp],[0,dir_fac*offset], linewidth=linewidth, color=color)
	l2 = Line2D([start_bp,start_bp+dir_fac*extent-dir_fac*(arrow_len*0.5)],[dir_fac*offset,dir_fac*offset], linewidth=linewidth, color=color)
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
	l1 = Line2D([start_bp,start_bp],[0,dir_fac*offset], linewidth=linewidth, color=color)
	l2 = Line2D([start_bp-extent,start_bp+extent],[dir_fac*offset,dir_fac*offset], linewidth=linewidth, color=color)
	ax.add_line(l1)
	ax.add_line(l2)
	# Shade the terminator area (normally smaller than symbol extent)
 	p2 = Polygon([(start_bp, -t_width), 
 		          (start_bp, t_width),
 		          (end_bp, t_width),
 		          (end_bp, -t_width)], facecolor=color, edgecolor=color, linewidth=linewidth, zorder=10)
	ax.add_patch(p2)

def draw_cds (ax, start_bp, end_bp, color=(0.0,1.0,0.0), linewidth=2, arrow_extent=16, body_offset=3, arrow_offset=2, 
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
		          facecolor=color, edgecolor=(0.0,0.0,0.0), linewidth=linewidth, hatch=hatch, zorder=10)
	ax.add_patch(p1)

def plot_library_arch (fig, gcl):
	return None

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
			draw_cds(ax, part_start_bp, part_end_bp, linewidth=linewidth, arrow_extent=100)
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

def plot_traces_with_arch (ax_arch, ax_traces, gcl, variant, traces, start_idx=None, end_idx=None):
	return None

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

# Example plot testing full single variant rendering
fig = plt.figure(figsize=(14,3))
ax = fig.add_subplot(3,1,1)
plot_variant_arch(ax, nifs, '25', start_idx=1, end_idx=-1, linewidth=1.2)

ax = fig.add_subplot(3,1,2)
plot_variant_arch(ax, nifs, '10', start_idx=1, end_idx=-1, linewidth=1.2)

ax = fig.add_subplot(3,1,3)
plot_variant_arch(ax, nifs, '1', start_idx=1, end_idx=-1, linewidth=1.2)

# Show the examples
plt.show()

