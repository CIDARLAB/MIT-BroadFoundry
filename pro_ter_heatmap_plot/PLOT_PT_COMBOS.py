from pt_heatmap_plotter import heatmap_plot
from load_values import load_value_array
from cluster_data import cluster_data

'''
This script creates a heatmap with promoters on the x-axis and terminators on the y.
'''

# This function reads in the data to be heatmapped
l_of_l = load_value_array('\Users\Eric\Dropbox (MIT)\Data\Flow',
                          'P-T Combo All Values',
                          'Log Sorted Vals')

data = l_of_l
#data = cluster_data(l_of_l)

heatmap_plot(data)
