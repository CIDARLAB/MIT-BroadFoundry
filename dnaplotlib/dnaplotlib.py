#!/usr/bin/env python
"""
dnaplotlib
==========
    This module is designed to allow for highly customisable visualisation of DNA
    fragments. Diagrams can be in the form of conceptual SBOL compliant icons or
    make use of icons whose width is scaled to allow for easier comparison of part 
    locations to trace information, such as for corresponding RNA-seq read depth
    data. All plotting is performed using matplotlib and to an axis object. This 
    enables the export of publication quality, vector-based figures. Furthermore,
    all standard renderers can be replaced with user defined versions to allow 
    for full customisation of the plot.
"""
#    dnaplotlib
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    Emerson Glassey <eglassey@mit.edu>
#    Bryan Der <bder@mit.edu>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

from matplotlib.patches import Polygon, Ellipse, Wedge, Circle
from matplotlib.lines   import Line2D
from math import sqrt

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT\n\
               Emerson Glassey <eglassey@mit.edu>, Voigt Lab, MIT\n\
               Bryan Der <bder@mit.edu>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

class DNARenderer:

	# Standard part types
	STD_PART_TYPES = ['Promoter',
	                  'CDS',
                      'Terminator',
                      'RBS',
                      'Scar',
                      'Spacer',
                      'Ribozyme',
                      'Ribonuclease',
                      'ProteinStability',
                      'Protease',
                      'Operator',
                      'Origin',
                      'Insulator']

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
			- part_renderers: standard renderer functions dict to use
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

					#update start,end for regulation
					part['start'] = prev_start
					part['end'] = prev_end

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
						
						#update start,end for regulation
						part['start'] = prev_start
						part['end'] = prev_end
						
						if first_part == True:
							first_start = prev_start
							first_part = False
			part_num += 1
		# Plot the regulatory links on the axis
		if regs != None:
			#return first_start, prev_end

			reg_num = 0
			for reg in regs:
				keys = reg.keys()

				# Check the part has minimal details required
				if 'type' in keys and 'start_part' in keys and 'end_part' in keys:
					# Extract custom part options (if available)

					reg_opts = None
					if 'opts' in reg.keys():
						reg_opts = reg['opts']
					
					if reg['type'] in reg_renderers.keys():
						reg_renderers[reg['type']](ax, reg['type'], 
							           reg_num, reg['start_part'], 
							           reg['end_part'], self.y_scale, 
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
	return stick_figure(ax,type,num,start,end,prev_end,y_scale,linewidth,opts)	
def sbol_protein_stability (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	return stick_figure(ax,type,num,start,end,prev_end,y_scale,linewidth,opts)	
def sbol_protease (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	return stick_figure(ax,type,num,start,end,prev_end,y_scale,linewidth,opts)
def sbol_ribonuclease (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	return stick_figure(ax,type,num,start,end,prev_end,y_scale,linewidth,opts)

def stick_figure (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 5.0
	y_extent = 10.0
	linestyle = '-'
	
	linetype  = "";
	shapetype = "";
	if(type == "Ribozyme"):
		linetype = 'dash'
		headgroup = 'O'
	elif(type == "Protease"):
		linetype = 'dash'
		headgroup = 'X'
	elif(type == "ProteinStability"):
		linetype = 'full'
		headgroup = 'O'
	elif(type == "Ribonuclease"):
		linetype = 'full'
		headgroup = 'X'

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
		x1 = Line2D([start,end],[-1*y_extent*1.25,-1*y_extent/1.5], 
		        	linewidth=linewidth, color=color, zorder=12, linestyle='-')
		x2 = Line2D([start,end],[-1*y_extent/1.5,-1*y_extent*1.25], 
		        	linewidth=linewidth, color=color, zorder=12, linestyle='-')

		dash1 = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,-y_extent/4], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		dash2 = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[-y_extent/2,-y_extent+(x_extent/2.0)], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		fullO = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,-y_extent+(x_extent/2.0)], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		fullX = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,-y_extent], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)

		if(headgroup == "O" and linetype == "dash"):
			ax.add_patch(c1)
			ax.add_line(dash1)
			ax.add_line(dash2)
		elif(headgroup == "X" and linetype == "dash"):
			ax.add_line(x1)
			ax.add_line(x2)
			ax.add_line(dash1)
			ax.add_line(dash2)
		elif(headgroup == "O" and linetype == "full"):
			ax.add_patch(c1)
			ax.add_line(fullO)
		elif(headgroup == "X" and linetype == "full"):
			ax.add_line(x1)
			ax.add_line(x2)
			ax.add_line(fullX)
		
	else:
		start = prev_end+start_pad
		end = start+x_extent
		final_end = end+end_pad
		rbs_center = (start+((end-start)/2.0),y_extent)
		c1 = Circle(rbs_center, x_extent/2.0, linewidth=linewidth, edgecolor=color, 
			        facecolor=(1,1,1), zorder=8)
		x1 = Line2D([start,end],[y_extent*1.25,y_extent/1.5], 
		        	linewidth=linewidth, color=color, zorder=12, linestyle='-')
		x2 = Line2D([start,end],[y_extent/1.5,y_extent*1.25], 
		        	linewidth=linewidth, color=color, zorder=12, linestyle='-')

		dash1 = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,y_extent/4], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		dash2 = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[y_extent/2,y_extent-(x_extent/2.0)], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		fullO = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,y_extent-(x_extent/2.0)], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)
		fullX = Line2D([end+((start-end)/2.0),end+((start-end)/2.0)],[0,-y_extent], 
			        linewidth=linewidth, color=color, zorder=8, linestyle=linestyle)

		if(headgroup == 'O' and linetype == 'dash'):
			ax.add_patch(c1)
			ax.add_line(dash1)
			ax.add_line(dash2)
		elif(headgroup == "X" and linetype == "dash"):
			ax.add_line(x1)
			ax.add_line(x2)
			ax.add_line(dash1)
			ax.add_line(dash2)
		elif(headgroup == "O" and linetype == "full"):
			ax.add_patch(c1)
			ax.add_line(fullO)
		elif(headgroup == "X" and linetype == "full"):
			ax.add_line(x1)
			ax.add_line(x2)
			ax.add_line(fullX)
		
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end


def sbol_scar (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 3.0
	y_extent = 1.0
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

	start = prev_end+start_pad
	end = start+x_extent
	final_end = end+end_pad
	
	l_top    = Line2D([start-x_extent,start+x_extent],[y_extent,y_extent], 
		        linewidth=linewidth, color=color, zorder=12, linestyle=linestyle)
	l_bottom = Line2D([start-x_extent,start+x_extent],[-1*y_extent,-1*y_extent], 
		        linewidth=linewidth, color=color, zorder=12, linestyle=linestyle)
	#white rectangle overlays backbone line
	p1 = Polygon([(start-x_extent, y_extent), 
		          (start-x_extent, -y_extent),
		          (start+x_extent, -y_extent),
		          (start+x_extent, y_extent)],
		          edgecolor=(1,1,1), facecolor=(1,1,1), linewidth=linewidth, zorder=11)		

	ax.add_patch(p1)
	ax.add_line(l_top)
	ax.add_line(l_bottom)

	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_spacer (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 6.0
	y_extent = 6.0
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
	
	start = prev_end+start_pad
	end = start+x_extent
	final_end = end+end_pad
	rbs_center = (start+((end-start)/2.0),0)
	center_x = start+(end-start)/2.0
	radius = x_extent/2

	delta = radius - 0.5 * radius * sqrt(2)

	l1 = Line2D([start+delta,end-delta],[radius-delta,-1*radius+delta], 
		        linewidth=linewidth, color=color, zorder=12, linestyle=linestyle)
	l2 = Line2D([start+delta,end-delta],[-1*radius+delta,radius-delta], 
		        linewidth=linewidth, color=color, zorder=12, linestyle=linestyle)
	c1 = Circle(rbs_center, x_extent/2.0, linewidth=linewidth, edgecolor=color, 
		        facecolor=(1,1,1), zorder=12)
	
	ax.add_patch(c1)
	ax.add_line(l1)
	ax.add_line(l2)

	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end


def sbol_origin (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 10.0
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
	
	start = prev_end+start_pad
	end = start+x_extent
	final_end = end+end_pad
	rbs_center = (start+((end-start)/2.0),0)
	
	c1 = Circle(rbs_center, x_extent/2.0, linewidth=linewidth, edgecolor=color, 
		        facecolor=(1,1,1), zorder=12)
	
	ax.add_patch(c1)
	
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_operator (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 3.0
	y_extent = 3.0
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

	start = prev_end+start_pad
	end = start+x_extent
	final_end = end+end_pad
	
	#white rectangle overlays backbone line
	p1 = Polygon([(start-x_extent, y_extent), 
		          (start-x_extent, -y_extent),
		          (start+x_extent, -y_extent),
		          (start+x_extent, y_extent)],
		          edgecolor=(0,0,0), facecolor=(1,1,1), linewidth=linewidth, zorder=11)		

	ax.add_patch(p1)
	
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

def sbol_insulator (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (0,0,0)
	start_pad = 2.0
	end_pad = 2.0
	x_extent = 3.0
	y_extent = 3.0
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

	start = prev_end+start_pad
	end = start+x_extent
	final_end = end+end_pad
	
	#white rectangle overlays backbone line
	p1 = Polygon([(start-x_extent, y_extent), 
		          (start-x_extent, -y_extent),
		          (start+x_extent, -y_extent),
		          (start+x_extent, y_extent)],
		          edgecolor=(0,0,0), facecolor=(1,1,1), linewidth=linewidth, zorder=12)		

	p2 = Polygon([((start-x_extent)-x_extent/2,  y_extent+x_extent/2), 
		          ((start-x_extent)-x_extent/2, -y_extent-x_extent/2),
		          ((start+x_extent)+x_extent/2, -y_extent-x_extent/2),
		          ((start+x_extent)+x_extent/2,  y_extent+x_extent/2)],
		          edgecolor=(0,0,0), facecolor=(1,1,1), linewidth=linewidth, zorder=11)		

	ax.add_patch(p1)
	ax.add_patch(p2)
	
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end


def temporary_repressor (ax, type, num, start, end, prev_end, y_scale, linewidth, opts):
	# Default options
	color = (1.0,0.0,0.0)
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
	
	e1center = (start+((end-start)/2.0),0)
	e2center = (start+((end-start)/2.0)+x_extent/3.75,0)

	e1 = Ellipse(e1center, y_extent/2, y_extent, edgecolor=(0,0,0), facecolor=color, 
				linewidth=linewidth, fill=True, zorder=12)
	e2 = Ellipse(e2center, y_extent/2, y_extent, edgecolor=(0,0,0), facecolor=color, 
				linewidth=linewidth, fill=True, zorder=11)

	ax.add_patch(e1)
	ax.add_patch(e2)

	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

###############################################################################
# Regulation renderers
###############################################################################

def repress (ax, type, num, start_part, end_part, y_scale, linewidth, opts):

	print start_part['start'],start_part['end']
	print end_part['start'],end_part['end']

	start = (start_part['start'] + start_part['end']) / 2
	end   = (end_part['start']   + end_part['end']) / 2

	line_away   = Line2D([start,start],[0,2], 
		        linewidth=linewidth, color=(0,0,0), zorder=12, linestyle='-')
	line_across = Line2D([start,end],[2,2], 
		        linewidth=linewidth, color=(0,0,0), zorder=12, linestyle='-')
	line_toward = Line2D([end,end],[2,0], 
		        linewidth=linewidth, color=(0,0,0), zorder=12, linestyle='-')
	ax.add_line(line_away)
	ax.add_line(line_across)
	ax.add_line(line_toward)

def induce (ax, type, num, start, end, y_scale, linewidth, opts):
	print 'call induce renderer'

#reg_renderers[reg['type']](ax, reg['type'], 
#								           reg_num, reg['start'], 
#								           reg['end'], self.y_scale, 
#								           self.linewidth, opts=reg_opts)


###############################################################################
# Trace Icon Renderers (icon width corrisponds to trace data)
###############################################################################

# TODO - Convert from GeneClusterLibrary visualisation package, and integrate dnaplotlib


