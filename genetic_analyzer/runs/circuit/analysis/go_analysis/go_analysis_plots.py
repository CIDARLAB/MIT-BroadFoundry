#!/usr/bin/env python


__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import numpy as np
import csv

import scipy
import scipy.stats as stats
import pylab as plt
import scipy.cluster.hierarchy as sch
import matplotlib.ticker
import matplotlib.ticker as ticker

from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

fmt_axis_outline_width = 2.0
fmt_line_width = 2.0

go_x_range = [-7.9, 7.9]


def load_go_data (de_filename):
	de_data = []
	# Open file and ignore headers
	de_reader = csv.reader(open(de_filename, 'rb'), delimiter='\t')
	header = next(de_reader, None)
	# Load all the data
	for row in de_reader:
		if len(row) > 6:
			go_term = row[2]
			p_corr = float(row[6])
			ratio_study = float(row[3].split('/')[0])/float(row[3].split('/')[1])
			ratio_background = float(row[4].split('/')[0])/float(row[4].split('/')[1])
			delta_ratio = np.log2(ratio_study/ratio_background)
			de_data.append([go_term, delta_ratio, p_corr])
	return de_data

###############################################################################
# PLOT THE FIGURES
###############################################################################

# Load the raw DE results
go_up_678 = load_go_data ('./results/broken_flask_vs_tube.up.go.txt')
go_down_678 = load_go_data ('./results/broken_flask_vs_tube.down.go.txt')

labs = [x[0] for x in go_up_678]
bar_vals = [x[1] for x in go_up_678]

bar_pvals = np.array([x[2] for x in go_up_678])

fig = plt.figure(figsize=(3.5,5.61))

ax = plt.subplot(1,2,2)
ax.barh(np.arange(len(labs))+0.25, bar_vals, height=0.5, facecolor=(0,0,0 ))
for axis in ['top','right', 'left']:
	ax.spines[axis].set_visible(False)
	ax.set_yticks([])
for axis in ['bottom']:
	ax.spines[axis].set_linewidth(fmt_axis_outline_width)
ax.xaxis.set_ticks_position("bottom")
plt.axvline(x=0, color=(0,0,0), linewidth=fmt_line_width)
plt.gca().invert_yaxis()
ax.set_xlim(go_x_range)
ax.set_xticklabels([])

ax = plt.subplot(1,2,1)
ax.barh(np.arange(len(labs))+0.25, np.log10(bar_pvals), height=0.5, facecolor=(0,0,0))
for axis in ['top','right', 'left']:
	ax.spines[axis].set_visible(False)
	ax.set_yticks([])
for axis in ['bottom']:
	ax.spines[axis].set_linewidth(fmt_axis_outline_width)
ax.xaxis.set_ticks_position("bottom")
plt.axvline(x=0, color=(0,0,0), linewidth=fmt_line_width)
plt.gca().invert_yaxis()
ax.set_xticklabels([])
ax.set_xlim([0.1,-4.5])

plt.tight_layout()
plt.subplots_adjust(hspace=.0)
fig.savefig('./results/go_up_regulated.pdf', transparent=True)
plt.close('all')

#####################################################################

y_lims = [0.0,2.0]

labs = [x[0] for x in go_down_678]
bar_vals = [x[1] for x in go_down_678]

bar_pvals = np.array([x[2] for x in go_down_678])

fig = plt.figure(figsize=(3.5,1.03))

ax = plt.subplot(1,2,2)
ax.barh(np.arange(len(labs))+0.25, bar_vals, height=0.5, facecolor=(0,0,0 ))
for axis in ['top','right', 'left']:
	ax.spines[axis].set_visible(False)
	ax.set_yticks([])
for axis in ['bottom']:
	ax.spines[axis].set_linewidth(fmt_axis_outline_width)
ax.xaxis.set_ticks_position("bottom")
plt.axvline(x=0, color=(0,0,0), linewidth=fmt_line_width)
plt.gca().invert_yaxis()
ax.set_xlim(go_x_range)
ax.set_xticklabels([])
ax.set_ylim(y_lims)

ax = plt.subplot(1,2,1)
ax.barh(np.arange(len(labs))+0.25, np.log10(bar_pvals), height=0.5, facecolor=(0,0,0))
for axis in ['top','right', 'left']:
	ax.spines[axis].set_visible(False)
	ax.set_yticks([])
for axis in ['bottom']:
	ax.spines[axis].set_linewidth(fmt_axis_outline_width)
ax.xaxis.set_ticks_position("bottom")
plt.axvline(x=0, color=(0,0,0), linewidth=fmt_line_width)
plt.gca().invert_yaxis()
ax.set_xticklabels([])
ax.set_xlim([0.1,-4.5])
ax.set_ylim(y_lims)

plt.tight_layout()
plt.subplots_adjust(hspace=.0)
fig.savefig('./results/go_down_regulated.pdf', transparent=True)
plt.close('all')


