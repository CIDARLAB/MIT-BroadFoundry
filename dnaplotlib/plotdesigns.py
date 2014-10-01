#!/usr/bin/env python
"""
	plotdeaigns.py

	Will plot the design of DNA constructs using SBOL notation.

 	Usage:
    ------
    python plotdesigns.py PARAM_FILENAME PART_FILENAME DESIGN_FILENAME OUT_FILENAME
"""
#    Plot Designs
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import sys
import getopt
import csv
import dnaplotlib as dpl
import matplotlib.pyplot as plt

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def make_float_if_needed (s):
	try:
		float(s)
		return float(s)
	except ValueError:
		return s

def load_part_information (filename):
	part_info = {}
	parts_reader = csv.reader(open(filename, 'rU'), delimiter=',')
	header = next(parts_reader)
	header_map = {}
	for i in range(len(header)):
		header_map[header[i]] = i
	attrib_keys = [k for k in header_map.keys() if k not in ['part_name', 'type']]
	for row in parts_reader:
		# Make the attributes map
		part_attribs_map = {}
		for k in attrib_keys:
			if row[header_map[k]] != '':
				if k == 'color':
					part_attribs_map[k] = [float(x) for x in row[header_map[k]].split(';')]
				else:
					part_attribs_map[k] = make_float_if_needed(row[header_map[k]])
		part_name = row[header_map['part_name']]
		part_type = row[header_map['type']]
		part_info[part_name] = [part_type, part_attribs_map]
	return part_info

def load_plot_parameters (filename):
	plot_params = {}
	param_reader = csv.reader(open(filename, 'rU'), delimiter=',')
	# Ignore header
	header = next(param_reader)
	# Process all parameters
	for row in param_reader:
		if len(row) >= 2:
			plot_params[row[0]] = make_float_if_needed(row[1])
	return plot_params

def load_dna_designs (filename, part_info):
	dna_designs = {}
	design_reader = csv.reader(open(filename, 'rU'), delimiter=',')
	# Ignore header
	header = next(design_reader)
	# Process all parameters
	for row in design_reader:
		if len(row[0]) != '':
			part_list = []
			for i in range(1,len(row)):
				# Handle reverse parts
				fwd = True
				part_name = row[i]
				if len(part_name) != 0:
					if part_name[0] == 'r':
						part_name = part_name[1:]
						fwd = False
					# Store the design
					part_design = {}
					cur_part_info = part_info[part_name]
					part_design['type'] = cur_part_info[0]
					if fwd == True:
						part_design['start'] = i
						part_design['end'] = i+1
					else:
						part_design['end'] = i
						part_design['start'] = i+1
					part_design['opts'] = cur_part_info[1]
					part_list.append(part_design)
			dna_designs[row[0]] = part_list
	return dna_designs

def plot_dna (dna_designs, out_filename, plot_params):
	# Create the renderer
	dr = dpl.DNARenderer(y_scale=plot_params['y_scale'], linewidth=plot_params['linewidth'])
	# We default to the SBOL part renderers
	part_renderers = {'Promoter'  :dpl.sbol_promoter, 
	                  'CDS'       :dpl.sbol_cds, 
                      'Terminator':dpl.sbol_terminator,
                      'RBS'       :dpl.sbol_rbs}
    # Create the figure
	fig = plt.figure(figsize=(plot_params['fig_y'],plot_params['fig_x']))
	# Cycle through the designs an plot on individual axes
	design_list = sorted(dna_designs.keys())
	num_of_designs = len(design_list)
	ax_list = []
	max_dna_len = 0.0
	for i in range(num_of_designs):
		# Create axis for the design and plot
		design =  dna_designs[design_list[i]]
		ax = fig.add_subplot(num_of_designs,1,i+1)
		if 'show_title' in plot_params.keys() and plot_params['show_title'] == 'Y':
			ax.set_title(design_list[i], fontsize=8)
		start, end = dr.renderDNA(ax, design, part_renderers)
		dna_len = end-start
		if max_dna_len < dna_len:
			max_dna_len = dna_len
		ax_list.append(ax)
	for ax in ax_list:
		ax.set_xticks([])
		ax.set_yticks([])
		# Set bounds
		ax.set_xlim([-0.01*max_dna_len,max_dna_len+(0.01*max_dna_len)])
		ax.set_ylim([-15,15])
		ax.set_aspect('equal')
		ax.set_axis_off()
	# Save the figure
	plt.tight_layout()
	fig.savefig(out_filename)
	# Clear the plotting cache
	plt.close('all')

def main():
	# parse command line options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)
	# process options
	for o, a in opts:
		if o in ("-h", "--help"):
			print __doc__
			sys.exit(0)
	# process arguments
	plot_params = load_plot_parameters(args[0])
	part_info = load_part_information(args[1])
	dna_designs = load_dna_designs (args[2], part_info)
	plot_dna(dna_designs, args[3], plot_params)

if __name__ == "__main__":
 	main()
