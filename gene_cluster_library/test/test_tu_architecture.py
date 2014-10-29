
import sys
sys.path.append('../')
import gene_cluster_library as gct

nifs = gct.GeneClusterLibrary()
nifs.load('./data/clean_nif_stata_library.txt')



# Extract the TUs

tus_basic = nifs.transcriptional_units(non_terminated=False, read_through=False)

tus_all = nifs.transcriptional_units(read_through=True)


