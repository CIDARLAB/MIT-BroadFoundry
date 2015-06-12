import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib import rcParams

rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['mathtext.fontset'] = 'stixsans'

from colorbar_settings import colorbar_legend_settings

'''
This script creates a heatmap with promoters on the x-axis and terminators on the y.
'''

def heatmap_plot(data) :

    array = data[0]
    x_title = 'Promoters'
    x_labels = data[1]
    y_title = 'Terminators'
    y_labels = data[2]
    
    fig, ax = plt.subplots()

    # These are mostly custom settings that are messy
    # and really cloud the formatting settings for the plot.
    # I moved them into a helper function for better readability
    # of the other code.
    cbarD = colorbar_legend_settings(5, 'log')

    im = plt.imshow(array,
                    interpolation='nearest',
                    cmap=cbarD['color_map'])


    # Format colorbar to show log10 demarcations and ticks outward, including minor tick marks
    cbar = fig.colorbar(im, ticks=cbarD['ticks'])
    
    cbar.ax.set_yticklabels(cbarD['scale'], fontsize=15)
    
    cbar.ax.tick_params(axis='y', direction='out')

    cbar.ax.yaxis.set_ticks(im.norm(cbarD['minorticks']), minor=True)

    # Set colorbar to same height as figure
    ax.axis('tight')
    
    # Move left and top spines outward by 0 points
    ax.spines['left'].set_position(('outward', 0))
    ax.spines['top'].set_position(('outward', 0))
    
    # Option to hide the right and bottom spines
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    
    # Only show ticks on the left and top spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('top')

    ax.set_xlabel(x_title, fontsize=18)
    ax.xaxis.set_label_position('top')
    ax.set_ylabel(y_title, fontsize=18)
    
    # Reset tick mark location from default 
    ax.locator_params(axis = 'x', nbins = 43)
    ax.locator_params(axis = 'y', nbins = 37)

    # Reset tick mark direction from default
    ax.tick_params(axis='x', direction='out')
    ax.tick_params(axis='y', direction='out')

    # Rename the labels using the lists created from the data
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='left')
    ax.set_yticklabels(y_labels)
    
    plt.show()
