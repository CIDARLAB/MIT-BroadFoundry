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
		# Plot the regulatory links to the axis
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


def sbol_promoter (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0.0,0.0,0.0)
	offset = 5.0
	extent = 20.0
	arrow_width = 1.0
	arrow_len = 6.0
	# Reset defaults if included
	if opts != None:
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
	if start > end:
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
	if start > end:
		return end, start
	else:
		return start, end

def sbol_cds (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	body_offset = 3
	arrow_extent = 16
	arrow_offset = 3
	hatch = ''
	color = (1,0,0)
	dir_fac = 1.0
	if start > end:
		dir_fac = -1.0
	p1 = Polygon([(start, body_offset), 
		          (start, -body_offset),
		          (end-dir_fac*arrow_extent, -body_offset),
		          (end-dir_fac*arrow_extent, -body_offset-arrow_offset),
		          (end, 0),
		          (end-dir_fac*arrow_extent, body_offset+arrow_offset),
		          (end-dir_fac*arrow_extent, body_offset)],
		          facecolor=color, edgecolor=(0.0,0.0,0.0), linewidth=linewidth, hatch=hatch, zorder=11)
	ax.add_patch(p1)
	if start > end:
		return end, start
	else:
		return start, end




#############################################################################################
# BASIC TESTING WILL BE REMOVED EVENTUALLY
#############################################################################################

# Create the render and DNA description
dr = DNARenderer(y_scale=1.0, linewidth=2.0)
parts = [{'type':'Promoter', 'start':10, 'end':30},
         {'type':'CDS', 'start':31, 'end':300},
         {'type':'Promoter', 'start':310, 'end':340}]
part_renderers = {'Promoter':sbol_promoter,
                  'CDS':sbol_cds}

# Set up the figure
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(1,1,1)

# Render the DNA
start, end = dr.renderDNA(ax, parts, part_renderers)
dna_len = end-start

# Set bounds
ax.set_xlim([start-(0.05*dna_len),end+(0.05*dna_len)])
ax.set_ylim([-30,30])
ax.set_axis_off()

# Save the figure
plt.tight_layout()
fig.savefig('test.pdf')

# Clear the plotting cache
plt.close('all')

