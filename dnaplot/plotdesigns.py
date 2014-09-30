#!/usr/bin/env python
"""
plotdesigns
===========
    Will plot designs of constructs held in CSV files
"""
#    Plot Designs
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import dnaplot as dp
import matplotlib.pyplot as plt

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def load_part_information (filename):
	part_info = {}



	return part_info

def load_plot_parameters (filename):

	plot_params = {}

	return plot_params

def load_dna_designs (filename, part_info):
	dna_designs = {}
	return dna_designs

def plot_dna (dna_designs, out_filename, plot_params):
	# Create the renderer
	dr = dp.DNARenderer(y_scale=plot_params['y_scale'], linewidth=plot_params['linewidth'])
	# We default to the SBOL part renderers
	part_renderers = {'Promoter'  :dp.sbol_promoter, 
	                  'CDS'       :dp.sbol_cds, 
                      'Terminator':dp.sbol_terminator,
                      'RBS'       :dp.sbol_rbs}
    # Create the figure
	fig = plt.figure(figsize=(plot_params['fig_x'],plot_params['fig_y']))
	# Cycle through the designs an plot on individual axes
	design_list = dna_designs.keys()
	num_of_designs = len(design_list)
	for i in range(num_of_designs):
		# Create axis for the design and plot
		design =  dna_designs[design_list[i]]
		ax = fig.add_subplot(1,num_of_designs,i)
		start, end = dr.renderDNA(ax, design, part_renderers)
		dna_len = end-start
		# Set bounds
		ax.set_xlim([start-(0.01*dna_len),end+(0.01*dna_len)])
		ax.set_ylim([-15,15])
		ax.set_aspect('equal')
		ax.set_axis_off()
	# Save the figure
	plt.tight_layout()
	fig.savefig(out_filename)
	# Clear the plotting cache
	plt.close('all')





parts = [{'type':'Promoter', 'start':0, 'end':30, 'opts':{'end_pad':-5}},
         {'type':'RBS', 'start':20, 'end':40},
         {'type':'CDS', 'start':250, 'end':300},
         {'type':'CDS', 'start':300, 'end':350, 'opts':{'hatch':'/////'}},
         {'type':'Terminator', 'start':351, 'end':352},
         {'type':'Promoter', 'start':360, 'end':380}]

