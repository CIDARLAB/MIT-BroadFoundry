#!/usr/bin/env python
"""
	Plot growth rates.
"""

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import csv

matplotlib.rcParams['lines.dash_joinstyle']  = 'miter'
matplotlib.rcParams['lines.dash_capstyle']   = 'butt'
matplotlib.rcParams['lines.solid_joinstyle'] = 'miter'
matplotlib.rcParams['lines.solid_capstyle']  = 'projecting'
# Make text editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype']          = 42

OUTPUT_PREFIX = './plots_general'

def plot_growth_rates (growth_data_tube, growth_data_flask, filename_out):
	fig = plt.figure(figsize=(2.2,1.5))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])

	width = 0.325
	bar_pad = 0.0
	ind = np.arange(len(growth_data_tube[0]))-width
	ax.bar(ind, growth_data_tube[0], width, yerr=growth_data_tube[1], linewidth=0.8, color=(0.5,0.5,0.5), edgecolor=(0,0,0),
		   error_kw=dict(ecolor=(0,0,0), lw=0.8, capsize=1.2, capthick=0.8, zorder=10))

	ind = np.arange(len(growth_data_flask[0]))+bar_pad
	ax.bar(ind, growth_data_flask[0], width, yerr=growth_data_flask[1], linewidth=0.8, color=(1,1,1), edgecolor=(0,0,0),
		   error_kw=dict(ecolor=(0,0,0), lw=0.8, capsize=1.2, capthick=0.8, zorder=10))


	ax.set_xlim([-0.7, 7.7])
	ax.set_ylim([0, 100])
	ax.tick_params(axis='y', which='major', labelsize=8, pad=2.5, length=2, width=0.5)
	ax.tick_params(axis='x', which='major', labelsize=8, pad=3, length=2, width=0.5)
	ax.set_xticklabels([], visible=False)
	for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(0.8)
			ax.spines[axis].set_linewidth(0.8)
			ax.spines[axis].set_linewidth(0.8)
	plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.15)
	plt.savefig(filename_out, transparent=True)
	plt.close('all')


growth_data_flask = [[44, 45.17004, 52.67905, 52.88924, 64.64209, 65.2663, 74.52865, 69.92028],
                     [1, 1.364975, 4.641541, 5.235608, 6.805679, 8.696901, 18.12313, 17.21088]]

growth_data_tube = [[44.28745, 44.70032, 50.66848, 49.73416, 52.78443, 54.38025, 59.83904, 59.89109],
                    [1.02123932, 1.014463848, 3.342831432, 3.240729338, 1.676952908, 2.487711914, 7.789541735, 8.430673703]]

plot_growth_rates(growth_data_tube, growth_data_flask, OUTPUT_PREFIX+'/growth_rates_tube.pdf')

