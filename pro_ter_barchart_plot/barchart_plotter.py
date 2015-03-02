import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import setp
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator

rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'in'
rcParams['mathtext.fontset'] = 'stixsans'
rcParams['ytick.major.size'] = 10
rcParams['ytick.major.width'] = 1
rcParams['ytick.minor.size'] = 7
rcParams['ytick.minor.width'] = 1

from scale_formatting import scale_formatting

def barchart_plot(data,
                  n_series,
                  n_points_in_series,
                  legend_labels,
                  legend_title,
                  color_scale,
                  xaxis_title,
                  xaxis_labels,
                  ylabel,
                  scale,
                  kind):
    
    # Instantiate the figure
    fig, ax = plt.subplots()

    # This is a way to set the index for x-axis spacing.
    # It is based on the number of points in the series
    ind = np.arange(n_points_in_series)

    # This sets the space between x-axis points
    margin = 0.1

    # This sets the width of the x-axis points
    width = (1.-2.*margin)/n_series

    # This makes the bars for each data point
    for i in range(0,n_series) :
        
        # This gives a different color value for each series
        colors = color_scale[0]+color_scale[1]*i/n_series

        # Command for bar creation with appropriate scaling for size of data
        ax.bar(ind+(width*i), data[i], width, color = cm.jet(colors))
    
    # It matters what order the formatting steps are done!

    # X-axis formatting commands
    ax.set_xlabel(xaxis_title, fontsize=20)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange((width*n_series)/2, n_points_in_series))
    ax.set_xticklabels(xaxis_labels, rotation=0, horizontalalignment='center')
    ax.tick_params(axis='x', labelsize=15)

    # Y-axis formatting commands
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_ylim(0,scale)
    ax.set_yticks(np.arange(scale+1))

    scale_settings = scale_formatting(scale, kind)
    
    ax.set_yticklabels(scale_settings[0], fontsize=15)
    ax.tick_params(axis='y', labelright=True, labelsize=15)
    
    minorLocator = FixedLocator(scale_settings[1])
    ax.yaxis.set_minor_locator(minorLocator)
    
    leg = ax.legend(legend_labels,
                    loc='upper center',
                    ncol = n_series,
                    title = legend_title,
                    fontsize =15)
    
    setp(leg.get_title(),fontsize='20')
    
    ax.margins(0.04, 0, tight=False)
    
    plt.show()



    
