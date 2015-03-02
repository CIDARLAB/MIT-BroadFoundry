import numpy as np
import scipy.cluster.hierarchy as sch
    
def cluster_data(data) :

    #Reformat list of lists into a numpy array
    full_array = np.array(data[0])
    
    # This is a numpy array, so columns and rows can be omitted at will
    array = full_array[:,:-2]

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
    y_labels = ['']+[y_labels[i] for i in idx1]
    array = array[idx1,:]
    idx2 = Z2['leaves']
    x_labels = ['']+[x_labels[i] for i in idx2]
    array = array[:,idx2]

    return [array, x_labels, y_labels]


if __name__ == "__maine__" :
    maine()
