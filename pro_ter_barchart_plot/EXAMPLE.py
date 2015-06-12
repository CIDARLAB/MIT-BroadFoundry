from pro_ter_data_selection import pro_ter_data_selection
from labeled_data_dict_builder import labeled_data_dict_builder
from barchart_plotter import barchart_plot

    
mode = 'Promoters'
#mode = 'Terminators'

parts = ['P13','P12', 'P23', 'P35', 'P30', 'P15']
#parts = ['T1','T4','T22']

# sort 'Y/N' and part_name to sort by
sort = ['Y', 'P15']

color_scale = [0.2, 0.4]
#color_scale = [0.67, 0.4]


                                    # Directory
dataD = labeled_data_dict_builder('C:\Python27\MIT-BroadFoundry\pro_ter_barchart_plot',
                                    # File
                                    'example_data',
                                    # Sheet name
                                    'Log Sorted Vals',
                                    # Y names column
                                    'A',
                                    # X names row
                                    '1',
                                    # first data cell
                                    ['C','3'])



chartD = pro_ter_data_selection(mode, parts, dataD, sort)

data = chartD['data']

barchart_plot(data,
                len(data),
                len(data[0]),
                chartD['legend_labels'],
                chartD['legend_title'],
                color_scale,
                chartD['xaxis_title'],
                chartD['xaxis_labels'],
                'GFP Fluorescence\n(corrected for autofluorescence)',
                5,
                'log')
