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

'''
This script creates a heatmap with promoters on the x-axis and terminators on the y.
The input values must be saved as log10 values in Excel in .xlsx format.
'''


def maine():
    pt_heatmap_plotter()

def pt_heatmap_plotter():

    # This function reads in the data to be heatmapped
    data = load_value_array()

    # These are mostly custom settings that are messy
    # and really cloud the formatting settings for the plot.
    # I moved them into a helper function for better readability
    # of the other code.
    plot_settings = plot_specific_settings()

    # This function plots the data. Many settings can be found here.
    heatmap_plot(data, plot_settings)

def load_value_array():

    wb_filename = 'P-T Combo All Values'

    sheetname = 'Sq Log Sorted Vals' 
    
    wb = load_workbook('\Users\Eric\Dropbox (MIT)\Data\Flow\%(f)s.xlsx' % {'f':wb_filename},
                       data_only=True)
    
    ws = wb.get_sheet_by_name(sheetname)

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
                    if 0 <= cell.value <= 10:
                        z_values.append(cell.value)
                        gfp_values.append(cell.value)
                    else :
                        z_values.append(0)
                        gfp_values.append(0)

        if row > 2:
                value_array.append(gfp_values)

    return [value_array, x_labels, y_labels]

def plot_specific_settings():

    # How many logs will the values span?
    # What scale are the numbers in Excel in?
    n_logs = 5
    scale = 'log'
    
    # This is a placeholder for more complex custom colormaps.
    # When using a standard colormap, enter it in the heatmap_plot
    color_map = 'jetwz'

    settingsD = {}
    
    # Define a new colormap so that '0' values are white or use a standard colormap
    jetwz = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]

    jetwz_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', jetwz, N=256)

    cmapD = {'jetwz': jetwz_map}

    settingsD['color_map'] = cmapD[color_map]

    # Generate scale and tick marks for color code bar
    stD = {4: [[r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'],[0,1,2,3,4]],
           5: [[r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'],[0,1,2,3,4,5]]}

    settingsD['scale'] = stD[n_logs][0]
    settingsD['ticks'] = stD[n_logs][1]
    
    # Set proper minor tick marks for scale being used
    minor_tickD = {'log':[0.301029996,0.477121255,0.602059991,0.698970004,0.77815125,0.84509804,0.903089987,0.954242509,
                          1.301029996,1.477121255,1.602059991,1.698970004,1.77815125,1.84509804,1.903089987,1.954242509,
                          2.301029996,2.477121255,2.602059991,2.698970004,2.77815125,2.84509804,2.903089987,2.954242509,
                          3.301029996,3.477121255,3.602059991,3.698970004,3.77815125,3.84509804,3.903089987,3.954242509],
                          #4.301029996,4.477121255,4.602059991,4.698970004],
                   'log_on_linear' :[2,3,4,5,6,7,8,9,1,
                                     20,30,40,50,60,70,80,90,
                                     200,300,400,500,600,700,800,900,
                                     2000,3000,4000,5000,6000,7000,8000,9000,
                                     20000,30000,40000,50000,60000,70000,80000,90000]}
    
    settingsD['minorticks'] = minor_tickD[scale]

    return settingsD

def heatmap_plot(data, plotD):
    
    fig, ax = plt.subplots()

    im = plt.imshow(data[0],
                    interpolation='nearest',
                    cmap=plotD['color_map'])

    # Format colorbar to show log10 demarcations and ticks outward, including minor tick marks
    cbar = fig.colorbar(im, ticks=plotD['ticks'])
    
    cbar.ax.set_yticklabels(plotD['scale'], fontsize=15)
    
    cbar.ax.tick_params(axis='y', direction='out')

    cbar.ax.yaxis.set_ticks(im.norm(plotD['minorticks']), minor=True)

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
    ax.set_xticklabels(data[1], rotation=45, horizontalalignment='left')
    ax.set_yticklabels(data[2])
    
    plt.show()


if __name__ == "__maine__" :
    maine()
