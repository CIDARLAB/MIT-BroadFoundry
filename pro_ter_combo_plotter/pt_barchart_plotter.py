import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from openpyxl import *
from openpyxl.cell import *
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

def maine():
    pt_barchart_plotter('T')

def pt_barchart_plotter(mode):

    promoters = [14,19,31]

    terminators = [4,10,21]

    titleD = {'P' : ['Promoter', 'Terminators'],
                 'T' : ['Terminators', 'Promoters']}
    
    color_scaleD = {'P' : [0.2, 0.4],
                    'T' : [0.67, 0.4]}

    # This logic retrieves data in either promoters vs. terminators
    # or terminators vs. promoters
    # Returns [[legend labels], [x-axis values], [data]]
    if titleD[mode][0] is 'Promoter' :
        pt_list = promoters_data(promoters)
    else :
        pt_list = terminators_data(terminators)

    data = pt_list[2]
    n_series = len(data)
    n_points_in_series = len(data[0])

    legend_labels = pt_list[0]

    x_axis_labels = pt_list[1]

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
        colors = color_scaleD[mode][0]+color_scaleD[mode][1]*i/n_series

        # Command for bar creation with appropriate scaling for size of data
        ax.bar(ind+(width*i), data[i], width, color = cm.jet(colors))
    
    # It matters what order the formatting steps are done!

    # X-axis formatting commands
    ax.set_xlabel(titleD[mode][1], fontsize=20)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange((width*n_series)/2, n_points_in_series))
    ax.set_xticklabels(x_axis_labels, rotation=0, horizontalalignment='center')
    ax.tick_params(axis='x', labelsize=15)

    # Y-axis formatting commands
    ax.set_ylabel('GFP Fluorescence\n(corrected for autofluorescence)', fontsize=15)
    ax.set_ylim(0,5)
    ax.set_yticks(np.arange(6))
    scale_5 = [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$']
    scale_6 = [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$']
    ax.set_yticklabels(scale_5, fontsize=15)
    ax.tick_params(axis='y', labelright=True, labelsize=15)
    
    minorLocator = FixedLocator(minor_tick_list('log'))
    ax.yaxis.set_minor_locator(minorLocator)
    
    leg = ax.legend(legend_labels,
                    loc='upper center',
                    ncol = n_series,
                    title = titleD[mode][0],
                    fontsize =15)
    
    setp(leg.get_title(),fontsize='20')
    
    ax.margins(0.04, 0, tight=False)
    
    plt.show()

def minor_tick_list(kind) :
    
    scale = 5
    
    if kind is 'log' :
        if scale is 5 :
            minor_ticks = [0.301029996,0.477121255,0.602059991,0.698970004,0.77815125,0.84509804,0.903089987,0.954242509,
                           1.301029996,1.477121255,1.602059991,1.698970004,1.77815125,1.84509804,1.903089987,1.954242509,
                           2.301029996,2.477121255,2.602059991,2.698970004,2.77815125,2.84509804,2.903089987,2.954242509,
                           3.301029996,3.477121255,3.602059991,3.698970004,3.77815125,3.84509804,3.903089987,3.954242509,
                           4.301029996,4.477121255,4.602059991,4.698970004,4.77815125,4.84509804,4.903089987,4.954242509]
        if scale is 6 :
            minor_ticks = [0.301029996,0.477121255,0.602059991,0.698970004,0.77815125,0.84509804,0.903089987,0.954242509,
                           1.301029996,1.477121255,1.602059991,1.698970004,1.77815125,1.84509804,1.903089987,1.954242509,
                           2.301029996,2.477121255,2.602059991,2.698970004,2.77815125,2.84509804,2.903089987,2.954242509,
                           3.301029996,3.477121255,3.602059991,3.698970004,3.77815125,3.84509804,3.903089987,3.954242509,
                           4.301029996,4.477121255,4.602059991,4.698970004,4.77815125,4.84509804,4.903089987,4.954242509,
                           5.301029996,5.477121255,5.602059991,5.698970004,5.77815125,5.84509804,5.903089987,5.954242509]

    if kind is 'log_on_linear' :
        if scale is 5 :
            minor_ticks = [2,3,4,5,6,7,8,9,1,
                           20,30,40,50,60,70,80,90,
                           200,300,400,500,600,700,800,900,
                           2000,3000,4000,5000,6000,7000,8000,9000,
                           20000,30000,40000,50000,60000,70000,80000,90000]

        if scale is 6 :
            minor_ticks = [2,3,4,5,6,7,8,9,1,
                           20,30,40,50,60,70,80,90,
                           200,300,400,500,600,700,800,900,
                           2000,3000,4000,5000,6000,7000,8000,9000,
                           20000,30000,40000,50000,60000,70000,80000,90000,
                           200000,300000,400000,500000,600000,700000,800000,900000]
    return minor_ticks

def promoters_data(promoters) :
    
    data = data_builder()
    
    pro_data = []
    pro_labels = []
    ter_labels = []

    for promoter in promoters :
        pro_labels.append(data[0][promoter-1])
        
    count = 0
    
    while count <= 25 :
        
        samp_list = []
        
        for promoter in promoters :
            
            samp_list.append(data[2][count][promoter-1])

        ok = True

        for x in samp_list :
            if x is 0 :
                ok = False

        if ok is True :
            pro_data.append(samp_list)
            
            ter_labels.append(data[1][count])
            
        count += 1

    p_count = 0
    transposed_data = []

    while p_count < len(promoters) :
        
        t_count = 0
        tr_data = []
        
        while t_count < len(ter_labels) :
            
            tr_data.append(pro_data[t_count][p_count])
            t_count += 1
            
        p_count += 1
        
        transposed_data.append(tr_data)
    
    return [pro_labels, ter_labels, transposed_data]


def terminators_data(terminators) :

    data = data_builder()
    
    ter_data = []
    ter_labels = []
    pro_labels = []
    
    for terminator in terminators :
        ter_labels.append(data[1][terminator-1])

    count = 0
    
    while count <= 35 :
        
        samp_list = []
        
        for terminator in terminators :
            
            samp_list.append(data[2][terminator-1][count])

        ok = True

        for x in samp_list :
            if x is 0 :
                ok = False

        if ok is True :
            ter_data.append(samp_list)
            
            pro_labels.append(data[0][count])
            
        count += 1

    t_count = 0
    transposed_data = []

    while t_count < len(terminators) :
        
        p_count = 0
        tr_data = []
        
        while p_count < len(pro_labels) :
            
            tr_data.append(ter_data[p_count][t_count])
            p_count += 1
            
        t_count += 1
        
        transposed_data.append(tr_data)

    return [ter_labels, pro_labels, transposed_data]

def data_builder():
    
    wb = load_workbook('\Users\Eric\Dropbox (MIT)\Data\Flow\Analysis #4.xlsx', data_only=True)
    
    ws = wb.get_sheet_by_name('Log AutoFL Adj Values by 1-4')

    cell_array = ws.rows

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
                if col > 2 :
                    x_labels.append(cell.value)

            if row > 2 :
                if col is 1 :
                    y_labels.append(cell.value)

            if row > 2:              
                if col > 2:
                    if 0 <= cell.value <= 10 :
                        gfp_values.append(cell.value)
                    else :
                        gfp_values.append(0)

        if row > 2:
                value_array.append(gfp_values)


    return [x_labels, y_labels, value_array]
    
if __name__ == "__maine__" :
    maine()
