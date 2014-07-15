
# How the library is structured (at present!)
import gene_cluster_library as gcl

# Load the Stata nif library data
nifs = gcl.GeneClusterLibrary()
nifs.load('./data/nif_stata_library.txt')

nifs.save_variant_genbank('2', 'output.gb')
