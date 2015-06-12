import numpy as np
import scipy.cluster.hierarchy as sch

def maine():
    from load_values import load_value_array

                            # Directory
    data = load_value_array('C:\Python27\MIT-BroadFoundry\pro_ter_heatmap_plot',
                            # File name
                            'example_data',
                            # Worksheet
                            'Log Sorted Vals')
    
    hardcode_cluster(data)
    #cluster_data(data)
    
def cluster_data(data) :

    #Reformat list of lists into a numpy array
    full_array = np.array(data[0])
    
    # This is a numpy array, so columns and rows can be omitted at will
    array = full_array[:,:]

    # Depending on format, adjust length of labels
    x_labels = data[1][1:]
    y_labels = data[2][1:]

    # This code does the dirty work of clustering
    Y = sch.linkage(array, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='right', labels=y_labels)
    Y2 = sch.linkage(np.transpose(array), method='centroid')
    Z2 = sch.dendrogram(Y2, orientation='top', labels=x_labels)

    # Sort D and genes on the clustering
    idx1 = Z1['leaves']
    print y_labels
    y_labels = ['']+[y_labels[i] for i in idx1]
    print y_labels
    array = array[idx1,:]
    idx2 = Z2['leaves']
    x_labels = ['']+[x_labels[i] for i in idx2]
    array = array[:,idx2]

    return [array, x_labels, y_labels]

def hardcode_cluster(data) :

    #Reformat list of lists into a numpy array
    full_array = np.array(data[0])
    
    # This is a numpy array, so columns and rows can be omitted at will
    array = full_array[:,:]

    # Depending on format, adjust length of labels
    x_labels = data[1][1:]
    y_labels = data[2][1:]


    idxx = []
    desired_x_ordering = ['P7','P8','P22','P24','P2','P32','P26','P16','P10',
                          'P23','P25','P5','P19','P9','P21','P18','P6','P20',
                          'P11','P1','P15','P30','P27','P12','P17','P28','P34',
                          'P4','P3','P14','P29','P33']
    
    for n, x in enumerate(desired_x_ordering):
        idxx.append(x_labels.index(x))
        
    x_labels = ['']+[x_labels[i] for i in idxx]
    array = array[:,idxx]
    
    idxy = []
    desired_y_ordering = ['T12','T14','T9','T7','T20','T19','T8','T11','T10','T16',
                          'T13','T5','T2','T18','T30','T6','T15','T17','T26','T22',
                          'T21','T1','T23','T4','T25','T24','T27','T3','T28','T29']

    for n, y in enumerate(desired_y_ordering):
        idxy.append(y_labels.index(y))

    y_labels = ['']+[y_labels[i] for i in idxy]
    array = array[idxy,:]


    return [array, x_labels, y_labels]


if __name__ == "__maine__" :
    maine()
