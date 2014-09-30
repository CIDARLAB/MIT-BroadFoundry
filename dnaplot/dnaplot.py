#!/usr/bin/env python
"""
dnaplot
=======
    This module is designed to allow for highly customisable visualisation of DNA
    fragments. Diagrams can be in the form of conceptual SBOL compliant icons or
    make use of scaling icons to allow for easier comparison of part locations to
    trace information, such as RNA-seq read depths. All plotting is performed using
    matplotlib to enable export of publication quality, vector-based figures. 
    Furthermore, all standard renderers can be replaced with user defined versions
    to allow for full customisation of the plot.
"""
#    DNA Plot
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

from matplotlib.patches import Polygon, Ellipse, Wedge
from matplotlib.lines import Line2D

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

class DNARenderer:

	# Standard part types
	STD_PART_TYPES = ['Promoter',
	                  'RBS',
	                  'CDS',
	                  'Terminator',
	                  'Scar',
	                  'Spacer',
	                  'Ribozyme']

	# Standard regulatory types
	STD_REG_TYPES = ['Repression',
	                 'Activation']

	def __init__(self, y_scale=1.0, linewidth=1.0):
		"""Constructor to generate an empty DNARenderer.
		"""
		self.y_scale = y_scale
		self.linewidth = linewidth

	def renderDNA(self, ax, parts, part_renderers, regs=None, reg_renderers=None):
		"""Render the parts on the DNA and regulation.
			- parts_list: list of dicts defining the parts
			- renderers: standard renderer functions dict to use
			- regs: list of regulations on the DNA
			- reg_renderers: dict of standard regulation renderers
		"""
		# Plot the parts to the axis
		part_num = 0
		prev_end = 0
		first_start = 0
		first_part = True
		for part in parts:
			keys = part.keys()
			# Check the part has minimal details required
			if 'type' in keys and 'start' in keys and 'end' in keys:
				# Extract custom part options (if available)
				part_opts = None
				if 'opts' in part.keys():
					part_opts = part['opts']
				# Use the correct renderer
				if 'renderer' in part.keys():
					# Use custom renderer
					prev_start, prev_end = part['renderer'](ax, part['type'], part_num, 
						             part['start'], part['end'], prev_end,
						             self.y_scale, self.linewidth, 
						             opts=part_opts)
					if first_part == True:
						first_start = prev_start
						first_part = False
				else:
					# Use standard renderer, if one exists
					if part['type'] in part_renderers.keys():
						prev_start, prev_end = part_renderers[part['type']](ax, 
							           part['type'], part_num, 
							           part['start'], part['end'], 
							           prev_end, self.y_scale, 
							           self.linewidth, opts=part_opts)
						if first_part == True:
							first_start = prev_start
							first_part = False
			part_num += 1
		# Plot the regulatory links on the axis
		if regs != None:
			reg_num = 0
			for reg in regs:
				keys = reg.keys()
				# Check the part has minimal details required
				if 'type' in keys and 'start' in keys and 'end' in keys:
					# Extract custom part options (if available)
					reg_opts = None
					if 'opts' in reg.keys():
						reg_opts = reg['opts']
					# Use the correct renderer
					if 'renderer' in part.keys():
						# Use custom renderer
						reg['renderer'](ax, reg['type'], reg_num, 
							            reg['start'], reg['end'], 
							            self.y_scale, self.linewidth,
							            opts=reg_opts)
					else:
						# Use standard renderer, if one exists
						if reg['type'] in reg_renderers.keys():
							reg_renderers[reg['type']](ax, reg['type'], 
								           reg_num, reg['start'], 
								           reg['end'], self.y_scale, 
								           self.linewidth, opts=reg_opts)
				reg_num += 1
		# Plot the backbone (z=1)
		l1 = Line2D([first_start,prev_end],[0,0], linewidth=self.linewidth, color=(0,0,0), zorder=1)
		ax.add_line(l1)
		return first_start, prev_end

###############################################################################
# SBOL Compliant Icons
###############################################################################

def sbol_promoter (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0.0,0.0,0.0)
	start_pad = 2.0
	end_pad = 2.0
	y_extent = 10
	x_extent = 10
	arrow_height = 2
	arrow_length = 4
	# Reset defaults if provided
	if opts != None:
		if 'color' in opts.keys():
			color = opts['color']
		if 'start_pad' in opts.keys():
			start_pad = opts['start_pad']
		if 'end_pad' in opts.keys():
			end_pad = opts['end_pad']
		if 'y_extent' in opts.keys():
			y_extent = opts['y_extent']
		if 'x_extent' in opts.keys():
			x_extent = opts['x_extent']
		if 'arrow_height' in opts.keys():
			arrow_height = opts['arrow_height']
		if 'arrow_length' in opts.keys():
			arrow_length = opts['arrow_length']
		if 'linewidth' in opts.keys():
			linewidth = opts['linewidth']
		if 'y_scale' in opts.keys():
			y_scale = opts['y_scale']
	# Check direction add start padding
	dir_fac = 1.0
	final_end = end
	final_start = prev_end
	if start > end:
		dir_fac = -1.0
		start = prev_end-start_pad
		end = start-x_extent
		final_end = end-end_pad
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
	# Draw the promoter symbol
	l1 = Line2D([start,start],[0,dir_fac*y_extent], linewidth=linewidth, color=color, zorder=9)
	l2 = Line2D([start,start+dir_fac*x_extent-dir_fac*(arrow_length*0.5)],[dir_fac*y_extent,dir_fac*x_extent], linewidth=linewidth, color=color, zorder=10)
	ax.add_line(l1)
	ax.add_line(l2)
	p1 = Polygon([(start+dir_fac*x_extent-dir_fac*arrow_length, dir_fac*y_extent+(arrow_height)), 
		          (start+dir_fac*x_extent, dir_fac*y_extent),
		          (start+dir_fac*x_extent-dir_fac*arrow_length, dir_fac*y_extent-(arrow_height))],
		          facecolor=color, edgecolor=color, linewidth=linewidth)
	ax.add_patch(p1)
	if final_start > final_end:
		return final_end, final_start
	else:
		return final_start, final_end

def sbol_cds (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (1,0,0)
	hatch = ''
	start_pad = 1.0
	end_pad = 1.0
	y_extent = 5
	x_extent = 50
	arrow_height = 4
	arrow_length = 8
	# Reset defaults if provided
	if opts != None:
		if 'color' in opts.keys():
			color = opts['color']
		if 'hatch' in opts.keys():
			hatch = opts['hatch']
		if 'start_pad' in opts.keys():
			start_pad = opts['start_pad']
		if 'end_pad' in opts.keys():
			end_pad = opts['end_pad']
		if 'y_extent' in opts.keys():
			y_extent = opts['y_extent']
		if 'x_extent' in opts.keys():
			x_extent = opts['x_extent']
		if 'arrow_height' in opts.keys():
			arrow_height = opts['arrow_height']
		if 'arrow_length' in opts.keys():
			arrow_length = opts['arrow_length']
		if 'linewidth' in opts.keys():
			linewidth = opts['linewidth']
		if 'y_scale' in opts.keys():
			y_scale = opts['y_scale']
	# Check direction add start padding
	dir_fac = 1.0
	final_end = end
	final_start = prev_end
	if start > end:
		dir_fac = -1.0
		start = prev_end-start_pad
		end = start-x_extent
		final_end = end-end_pad
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
	# Draw the CDS symbol
	p1 = Polygon([(start, y_extent), 
		          (start, -y_extent),
		          (end-dir_fac*arrow_length, -y_extent),
		          (end-dir_fac*arrow_length, -y_extent-arrow_height),
		          (end, 0),
		          (end-dir_fac*arrow_length, y_extent+arrow_height),
		          (end-dir_fac*arrow_length, y_extent)],
		          edgecolor=(0.0,0.0,0.0), facecolor=color, linewidth=linewidth, hatch=hatch, zorder=11)
	ax.add_patch(p1)
	if final_start > final_end:
		return final_end, final_start
	else:
		return final_start, final_end

def sbol_terminator (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	y_extent = 10.0
	x_extent = 4.0
	# Reset defaults if provided
	if opts != None:
		if 'color' in opts.keys():
			color = opts['color']
		if 'start_pad' in opts.keys():
			start_pad = opts['start_pad']
		if 'end_pad' in opts.keys():
			end_pad = opts['end_pad']
		if 'y_extent' in opts.keys():
			y_extent = opts['y_extent']
		if 'x_extent' in opts.keys():
			x_extent = opts['x_extent']
		if 'linewidth' in opts.keys():
			linewidth = opts['linewidth']
		if 'y_scale' in opts.keys():
			y_scale = opts['y_scale']
	# Check direction add start padding
	dir_fac = 1.0
	final_end = end
	final_start = prev_end
	if start > end:
		dir_fac = -1.0
		start = prev_end-start_pad
		end = start-x_extent
		final_end = end-end_pad
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
	# Draw the terminator symbol
	l1 = Line2D([start,start],[0,dir_fac*y_extent], linewidth=linewidth, color=color, zorder=8)
	l2 = Line2D([start-x_extent,start+x_extent],[dir_fac*y_extent,dir_fac*y_extent], linewidth=linewidth, color=color, zorder=9)
	ax.add_line(l1)
	ax.add_line(l2)
	if final_start > final_end:
		return final_end, final_start
	else:
		return final_start, final_end

def sbol_rbs (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,1)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 10.0
	# Reset defaults if provided
	if opts != None:
		if 'color' in opts.keys():
			color = opts['color']
		if 'start_pad' in opts.keys():
			start_pad = opts['start_pad']
		if 'end_pad' in opts.keys():
			end_pad = opts['end_pad']
		if 'x_extent' in opts.keys():
			x_extent = opts['x_extent']
		if 'linewidth' in opts.keys():
			linewidth = opts['linewidth']
		if 'y_scale' in opts.keys():
			y_scale = opts['y_scale']
	# Check direction add start padding
	dir_fac = 1.0
	final_end = end
	final_start = prev_end
	rbs_center = (0,0)
	if start > end:
		start = prev_end+start_pad
		end = prev_end+start_pad+x_extent
		final_end = end+end_pad
		rbs_center = (end+((start-end)/2.0),0)
		w1 = Wedge(rbs_center, x_extent/2.0, 180, 360, linewidth=linewidth, facecolor=color, zorder=8)
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
		rbs_center = (start+((end-start)/2.0),0)
		w1 = Wedge(rbs_center, x_extent/2.0, 0, 180, linewidth=linewidth, facecolor=color, zorder=8)
	# Draw the RBS symbol
	ax.add_patch(w1)
	if final_start > final_end:
		return final_end, final_start
	else:
		return final_start, final_end



#############################################################################################
# BASIC TESTING WILL BE REMOVED EVENTUALLY
#############################################################################################

# Create the renderer and DNA description
dr = DNARenderer(y_scale=1.0, linewidth=1.2)

parts = [{'type':'Promoter', 'start':0, 'end':30, 'opts':{'end_pad':-5}},
         {'type':'RBS', 'start':20, 'end':40},
         {'type':'CDS', 'start':250, 'end':300},
         {'type':'CDS', 'start':300, 'end':350, 'opts':{'hatch':'/////'}},
         {'type':'Terminator', 'start':351, 'end':352},
         {'type':'Promoter', 'start':360, 'end':380}]

parts_rev = [{'type':'Terminator', 'start':100, 'end':150},
		{'type':'CDS', 'start':150, 'end':200, 'opts':{'hatch':'/////'}},
         {'type':'RBS', 'start':200, 'end':250},
         {'type':'Promoter', 'start':250, 'end':300, 'opts':{'end_pad':-5}}]

part_renderers = {'Promoter':sbol_promoter, 'CDS':sbol_cds, 
                  'Terminator':sbol_terminator, 'RBS':sbol_rbs}

# Set up the figure
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(3,2))
ax = fig.add_subplot(1,1,1)
# Render the DNA
start, end = dr.renderDNA(ax, parts, part_renderers)
dna_len = end-start
# Set bounds
ax.set_xlim([start-(0.01*dna_len),end+(0.01*dna_len)])
ax.set_ylim([-15,15])
ax.set_aspect('equal')
ax.set_axis_off()
# Save the figure
plt.tight_layout()
fig.savefig('test_fwd.pdf')


fig = plt.figure(figsize=(3,2))
ax = fig.add_subplot(1,1,1)
# Render the DNA
start, end = dr.renderDNA(ax, parts, part_renderers)
dna_len = end-start
# Set bounds
ax.set_xlim([start-(0.01*dna_len),end+(0.01*dna_len)])
ax.set_ylim([-15,15])
ax.set_aspect('equal')
ax.set_axis_off()
# Save the figure
plt.tight_layout()
fig.savefig('test_rev.pdf')

# Clear the plotting cache
plt.close('all')

