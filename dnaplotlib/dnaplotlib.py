#!/usr/bin/env python
"""
dnaplotlib
==========
    This module is designed to allow for highly customisable visualisation of DNA
    fragments. Diagrams can be in the form of conceptual SBOL compliant icons or
    make use of icons whose wide is scaled to allow for easier comparison of part 
    locations to trace information, such as for corresponding RNA-seq read depth
    data. All plotting is performed using matplotlib and to an axis object. This 
    enables the export of publication quality, vector-based figures. Furthermore,
    all standard renderers can be replaced with user defined versions to allow 
    for full customisation of the plot.
"""
#    dnaplotlib
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

from matplotlib.patches import Polygon, Ellipse, Wedge, Circle
from matplotlib.lines   import Line2D

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

	def __init__(self, y_scale=1.0, linewidth=1.0, 
		         backbone_pad_left=0.0, backbone_pad_right=0.0):
		"""Constructor to generate an empty DNARenderer.
		"""
		self.y_scale = y_scale
		self.linewidth = linewidth
		self.backbone_pad_left = backbone_pad_left
		self.backbone_pad_right = backbone_pad_right

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
		l1 = Line2D([first_start-self.backbone_pad_left,prev_end+self.backbone_pad_right],[0,0], 
			        linewidth=self.linewidth, color=(0,0,0), zorder=10)
		ax.add_line(l1)
		return first_start, prev_end

###############################################################################
# SBOL Compliant Icon Renderers
###############################################################################

def sbol_promoter (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0.0,0.0,0.0)
	start_pad = 2.0
	end_pad = 2.0
	y_extent = 10
	x_extent = 10
	arrowhead_height = 2
	arrowhead_length = 4
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
		if 'arrowhead_height' in opts.keys():
			arrowhead_height = opts['arrowhead_height']
		if 'arrowhead_length' in opts.keys():
			arrowhead_length = opts['arrowhead_length']
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
		start = prev_end+end_pad+x_extent
		end = prev_end+end_pad
		final_end = start+start_pad
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
	# Draw the promoter symbol
	l1 = Line2D([start,start],[0,dir_fac*y_extent], linewidth=linewidth, 
		        color=color, zorder=9)
	l2 = Line2D([start,start+dir_fac*x_extent-dir_fac*(arrowhead_length*0.5)],
                [dir_fac*y_extent,dir_fac*x_extent], linewidth=linewidth, 
                color=color, zorder=10)
	ax.add_line(l1)
	ax.add_line(l2)
	p1 = Polygon([(start+dir_fac*x_extent-dir_fac*arrowhead_length, 
		           dir_fac*y_extent+(arrowhead_height)), 
		          (start+dir_fac*x_extent, dir_fac*y_extent),
		          (start+dir_fac*x_extent-dir_fac*arrowhead_length, 
		           dir_fac*y_extent-(arrowhead_height))],
		          facecolor=color, edgecolor=color, linewidth=linewidth)
	ax.add_patch(p1)
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_cds (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0.7,0.7,0.7)
	hatch = ''
	start_pad = 1.0
	end_pad = 1.0
	y_extent = 5
	x_extent = 30
	arrowhead_height = 4
	arrowhead_length = 8
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
		if 'arrowhead_height' in opts.keys():
			arrowhead_height = opts['arrowhead_height']
		if 'arrowhead_length' in opts.keys():
			arrowhead_length = opts['arrowhead_length']
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
		start = prev_end+end_pad+x_extent
		end = prev_end+end_pad
		final_end = start+start_pad
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
	# Draw the CDS symbol
	p1 = Polygon([(start, y_extent), 
		          (start, -y_extent),
		          (end-dir_fac*arrowhead_length, -y_extent),
		          (end-dir_fac*arrowhead_length, -y_extent-arrowhead_height),
		          (end, 0),
		          (end-dir_fac*arrowhead_length, y_extent+arrowhead_height),
		          (end-dir_fac*arrowhead_length, y_extent)],
		          edgecolor=(0.0,0.0,0.0), facecolor=color, linewidth=linewidth, 
		          hatch=hatch, zorder=11)
	ax.add_patch(p1)
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

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
		start = prev_end+end_pad+x_extent
		end = prev_end+end_pad
		final_end = start+start_pad
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
	# Draw the terminator symbol
	l1 = Line2D([start,start],[0,dir_fac*y_extent], linewidth=linewidth, 
		        color=color, zorder=8)
	l2 = Line2D([start-x_extent,start+x_extent],[dir_fac*y_extent,dir_fac*y_extent], 
		        linewidth=linewidth, color=color, zorder=9)
	ax.add_line(l1)
	ax.add_line(l2)
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_rbs (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
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
		start = prev_end+end_pad+x_extent
		end = prev_end+end_pad
		final_end = start+start_pad
		rbs_center = (end+((start-end)/2.0),0)
		w1 = Wedge(rbs_center, x_extent/2.0, 180, 360, linewidth=linewidth, 
			       facecolor=color, zorder=8)
		ax.add_patch(w1)
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
		rbs_center = (start+((end-start)/2.0),0)
		w1 = Wedge(rbs_center, x_extent/2.0, 0, 180, linewidth=linewidth, 
			       facecolor=color, zorder=8)
		ax.add_patch(w1)
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_ribozyme (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 5.0
	y_extent = 10.0
	linestyle = '-'
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
		if 'y_extent' in opts.keys():
			y_extent = opts['y_extent']
		if 'linestyle' in opts.keys():
			linestyle = opts['linestyle']
		if 'linewidth' in opts.keys():
			linewidth = opts['linewidth']
		if 'y_scale' in opts.keys():
			y_scale = opts['y_scale']
	# Check direction add start padding
	final_end = end
	final_start = prev_end
	if start > end:
		start = prev_end+end_pad+x_extent
		end = prev_end+end_pad
		final_end = start+start_pad
		rbs_center = (end+((start-end)/2.0),-y_extent)
		c1 = Circle(rbs_center, x_extent/2.0, linewidth=linewidth, edgecolor=color, 
			        facecolor=(1,1,1), zorder=8)
		ax.add_patch(c1)
		l1 = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,-y_extent+(x_extent/2.0)], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		ax.add_line(l1)
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
		rbs_center = (start+((end-start)/2.0),y_extent)
		c1 = Circle(rbs_center, x_extent/2.0, linewidth=linewidth, edgecolor=color, 
			        facecolor=(1,1,1), zorder=8)
		ax.add_patch(c1)
		l1 = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,y_extent-(x_extent/2.0)], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		ax.add_line(l1)
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

###############################################################################
# Trace Icon Renderers (icon width corrisponds to trace data)
###############################################################################

# TODO - Convert from GeneClusterLibrary visualisation package, and integrate dnaplotlib


