from pro_ter_data_selection import pro_ter_data_selection
from labeled_data_dict_builder import labeled_data_dict_builder
from barchart_plotter import barchart_plot

def maine():
    pt_barchart_plotter()

def pt_barchart_plotter():
    
    mode = 'Promoters'

    #mode = 'Terminators'

    parts = ['P5', 'P4']
    
    #parts = ['T1','T4','T22']

    color_scale = [0.2, 0.4]
    #color_scale = [0.67, 0.4]


                                      # Directory
    dataD = labeled_data_dict_builder('\Users\Eric\Dropbox (MIT)\Data\Flow',
                                      # File
                                      'P-T Combo All Values',
                                      # Sheet name
                                      'Log Sorted Vals',
                                      # Y names column
                                      'A',
                                      # X names row
                                      '1',
                                      # first data cell
                                      ['C','3'])



    chartD = pro_ter_data_selection(mode, parts, dataD)

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


    
if __name__ == "__maine__" :
    maine()
