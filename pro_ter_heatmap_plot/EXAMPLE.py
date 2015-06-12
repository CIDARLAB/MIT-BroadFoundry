from pt_heatmap_plotter import heatmap_plot
from load_values import load_value_array
from cluster_data import cluster_data

'''
This script creates a heatmap with promoters on the x-axis and terminators on the y.
'''

# This function reads in the data to be heatmapped
# The function returns a list of lists of the data

                                  # Directory
list_of_lists = load_value_array('C:\Python27\MIT-BroadFoundry\pro_ter_heatmap_plot',
                                  # File name
                                 'example_data',
                                  # Worksheet
                                 'Log Sorted Vals')


# Use list of list for data already sorted in Excel
data = list_of_lists

# Use cluster_data() for clustering using T. Gorochowski's
# Automated clustering functions - see cluster_data.py
#data = cluster_data(list_of_lists)

heatmap_plot(data)
