#!/usr/bin/env python
"""
dnaplot
=======

    This module is designed to allow for highly customisable visualisation of dna
    fragments. Diagrams can be in the form of conceptual SBOL compliant icons or
    make use of scaling icons to allow for easier comparison of part locations to
    trace information, such as RNA-seq read depths. All plotting is performed using
    matplotlib to enable export of publication quality, vector-based figures. 
    Furthermore, all standard renderers can be replaced with user defined versions
    to allow for full customisation of the plot.
"""
#    DNA plot
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

from matplotlib.patches import Polygon, Ellipse
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

	def __init__(self, y_scale=10, linewidth=1):
		"""Constructor to generate an empty DNARenderer.
		"""
		self.y_scale = y_scale
		self.linewidth = linewidth

	def renderDNA(ax, parts, part_renderers, regs=None, reg_renderers=None):
		"""Render the parts on the DNA and regulation.
			- parts_list: list of dicts defining the parts
			- renderers: standard renderer functions dict to use
			- regs: list of regulations on the DNA
			- reg_renderers: dict of standard regulation renderers
		"""
		# Plot the backbone (z=1)
		start = parts[0]['start']
		end = parts[-1]['end']

		l1 = Line2D([start,start],[0,dir_fac*offset], linewidth=self.linewidth, color=color, zorder=1)
		ax.add_line(l1)
		# Plot the parts to the axis
		part_num = 0
		for part in parts:
			keys = part.keys()
			# Check the part has minimal details required
			if 'type' in keys and 'start' in keys and 'end' in keys and 'direction' in keys:
				# Extract custom part options (if available)
				part_opts = None
				if 'opts' in part.keys():
					part_opts = part['opts']
				# Use the correct renderer
				if 'renderer' in part.keys():
					# Use custom renderer
					part['renderer'](ax, part['type'], part_num, 
						             part['start'], part['end'],
						             part['direction'], self.y_scale, 
						             self.linewidth, opts=part_opts)
				else:
					# Use standard renderer, if one exists
					if part['type'] in part_renderers.keys():
						part_renderers[part['type']](ax, 
							           part['type'], part_num, 
							           part['start'], part['end'], 
							           part['direction'], self.y_scale, 
							           self.linewidth, opts=part_opts)
			part_num += 1
		# Plot the regulatory links to the axis
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
		return ax


def sbol_promoter (ax, type, num, start, end, direction, y_scale, linewidth, opts):
	# Default options
	color = (1.0,0.0,0.0)
	offset = 1.0
	extent = 10.0
	arrow_width = 10.0
	arrow_len = 2.0
	# Reset defaults if included
	if 'color' in opts.keys():
		color = opts['color']
	if 'offset' in opts.keys():
		offset = opts['offset']
	if 'extent' in opts.keys():
		extent = opts['extent']
	if 'arrow_width' in opts.keys():
		arrow_width = opts['arrow_width']
	if 'arrow_len' in opts.keys():
		arrow_len = opts['arrow_len']
	if 'linewidth' in opts.keys():
		linewidth = opts['linewidth']
	if 'y_scale' in opts.keys():
		y_scale = opts['y_scale']
	# Draw the promoter symbol
	dir_fac = 1.0
	if direction = 'F':
		dir_fac = -1.0
	l1 = Line2D([start,start],[0,dir_fac*offset], linewidth=linewidth, color=color, zorder=9)
	l2 = Line2D([start,start+dir_fac*extent-dir_fac*(arrow_len*0.5)],[dir_fac*offset,dir_fac*offset], linewidth=linewidth, color=color, zorder=10)
	ax.add_line(l1)
	ax.add_line(l2)
	p1 = Polygon([(start+dir_fac*extent-dir_fac*arrow_len, dir_fac*offset+(arrow_width)), 
		          (start+dir_fac*extent, dir_fac*offset),
		          (start+dir_fac*extent-dir_fac*arrow_len, dir_fac*offset-(arrow_width))],
		          facecolor=color, edgecolor=color, linewidth=linewidth)
	ax.add_patch(p1)



# Create the render and DNA description
DNARenderer dr(y_scale=11, linewidth=2.0)
parts = [{'type':'Promoter', 'start':10, 'end':30, 'direction':'F'},
         {'type':'CDS', 'start':31, 'end':300, 'direction':'F'},
         {'type':'Promoter', 'start':310, 'end':340, 'direction':'F'}]
part_renderers = {'Promoter':sbol_promoter}

# Set up the figure



# Render the DNA
dr.renderDNA(ax, parts, part_renderers)





