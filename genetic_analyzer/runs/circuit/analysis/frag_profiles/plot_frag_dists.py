#!/usr/bin/env python
"""
Plot the fragment size distributions
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI OSL 3.0'
__version__ = '1.0'

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv

matplotlib.rcParams['lines.dash_joinstyle']  = 'miter'
matplotlib.rcParams['lines.dash_capstyle']   = 'butt'
matplotlib.rcParams['lines.solid_joinstyle'] = 'miter'
matplotlib.rcParams['lines.solid_capstyle']  = 'projecting'
# Make text editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype']          = 42 

def load_frag_file (filename):
	data = [[],[]]
	f = open(filename, 'rU')
	reader = csv.reader(f, delimiter='\t')
	for row in reader:
		if len(row) == 2:
			data[0].append(float(row[0]))
			data[1].append(float(row[1]))
	return np.array(data)

def plot_frag (ax, filename, show_y=True, show_x=True):
	shaded_len = 78
	data = load_frag_file(filename)
	ax.fill_between(data[0,:], data[1,:]/float(np.sum(data[1,:])), np.zeros(len(data[1,:])), linewidth=0.0, facecolor=(0.5,0.5,0.5))
	#plt.axvline(x=1, linewidth=0.8, color=(0,0,0))
	ax.fill_between([0,shaded_len], [0,0], [1,1], linewidth=0.0, facecolor=(0.7,0.7,0.7), alpha=0.5)
	ax.tick_params(axis='x', labelsize=10)
	ax.tick_params(axis='y', labelsize=10)
	if show_y == False:
		ax.set_yticklabels([], visible=False)
	else:
		ax.set_yticklabels(['0.000', '', '0.002', '', '0.004', '', '0.006'])
		ax.set_ylabel('Density', fontsize=10, labelpad=2)
	if show_x == False:
		ax.set_xticklabels([], visible=False)
	else:
		ax.set_xticklabels(['',100, '', 300, '', 500])
		ax.set_xlabel('Length (bp)', fontsize=10, labelpad=1)
	for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(1.0)
	ax.set_xlim([0,600])
	ax.set_ylim([0,0.007])
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.tick_params('both', length=2.5, width=0.5, which='major')
	
def plot_frag_dist (cond, out_filename):
	gs = gridspec.GridSpec(2, 4)
	fig = plt.figure(figsize=(6,2.5))
	labs = ['-/-/-', '+/-/-', '-/+/-', '+/+/-', 
        '-/-/+', '+/-/+', '-/+/+', '+/+/+']
	for i in range(8):
		ax = plt.subplot(gs[i])
		show_y = True
		show_x = True
		if i not in [0,4]:
			show_y = False
		if i in range(4):
			show_x = False
		plot_frag (ax, './data/Circuit_sample_T'+str(i+1)+'_'+cond+'_span_dist.txt', show_x=show_x, show_y=show_y)
		ax.text(0.92,0.92,labs[i], horizontalalignment='right', verticalalignment='top',transform=ax.transAxes, fontsize=10)
	plt.subplots_adjust(left=0.1, right=0.98, top=0.94, bottom=0.15, wspace=0.05, hspace=0.08)
	plt.savefig(out_filename, transparent=True)


plot_frag_dist ('1', 'frag_dist_1_tube.pdf')
plot_frag_dist ('2', 'frag_dist_2_flask.pdf')

