import numpy as np
from openpyxl import *
from openpyxl.cell import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FixedLocator
from matplotlib import rcParams

rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['mathtext.fontset'] = 'stixsans'


def maine():
    pt_heatmap_plotter()

def pt_heatmap_plotter():
    
    wb = load_workbook('\Users\Eric\Dropbox (MIT)\Data\Flow\Analysis #4.xlsx', data_only=True)
    
    ws = wb.get_sheet_by_name('Log New Values')

    cell_array = ws.rows

    z_values = []
    x_labels = []
    y_labels = []

    value_array = []

    for row in cell_array :
        gfp_values = []
        for cell in row :
            xy = coordinate_from_string(cell.coordinate)
            col = column_index_from_string(xy[0])
            row = xy[1]

            if row is 1 :
                if col >= 2 :
                    x_labels.append(cell.value)

            if row >= 2 :
                if col is 1 :
                    y_labels.append(cell.value)

            if row > 2:              
                if col > 2:
                    if 0 <= cell.value <= 10 :
                        z_values.append(cell.value)
                        gfp_values.append(cell.value)
                    else :
                        z_values.append(0)
                        gfp_values.append(0)

        if row > 2:
                value_array.append(gfp_values)
   
    fig, ax = plt.subplots()

    # Define a new colormap so that '0' values are white
    colors = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]
    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

    im = plt.imshow(value_array, interpolation='nearest', cmap=new_map)

    # Format colorbar to show log10 demarcations and ticks outward
    cbar = fig.colorbar(im, ticks=[0,1,2,3,4,5])
    cbar.ax.set_yticklabels([r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'], fontsize=15)
    cbar.ax.tick_params(axis='y', direction='out')

    minorticks = im.norm(minor_tick_list('log'))
    cbar.ax.yaxis.set_ticks(minorticks, minor=True)

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

    ax.set_xlabel('Promoters', fontsize=18)
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Terminators', fontsize=18)
    
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

def minor_tick_list(kind) :

    if kind is 'log' :
        minor_ticks = [0.301029996,0.477121255,0.602059991,0.698970004,0.77815125,0.84509804,0.903089987,0.954242509,
                       1.301029996,1.477121255,1.602059991,1.698970004,1.77815125,1.84509804,1.903089987,1.954242509,
                       2.301029996,2.477121255,2.602059991,2.698970004,2.77815125,2.84509804,2.903089987,2.954242509,
                       3.301029996,3.477121255,3.602059991,3.698970004,3.77815125,3.84509804,3.903089987,3.954242509,
                       4.301029996,4.477121255,4.602059991,4.698970004]

    if kind is 'log_on_linear' :
        minor_ticks = [2,3,4,5,6,7,8,9,1,
                       20,30,40,50,60,70,80,90,
                       200,300,400,500,600,700,800,900,
                       2000,3000,4000,5000,6000,7000,8000,9000,
                       20000,30000,40000,50000,60000,70000,80000,90000]

    return minor_ticks

    
if __name__ == "__maine__" :
    maine()
