
import sys
sys.path.append('../')
import gene_cluster_library as gct
import gene_cluster_analysis as gca

# Load the example data set
nifs = gct.GeneClusterLibrary()
nifs.load('./data/clean_nif_stata_library.txt')

# Extract the TUs (test functions don't error)
tus_basic = nifs.transcriptional_units(non_terminated=False, read_through=False)
tus_all = nifs.transcriptional_units(read_through=True)
print 'There are', len(tus_basic), 'normal TUs and', len(tus_all), 'full.'

# Extract all the sequences for the TUs (we can use this data for fitting)
gca.generate_all_variant_tus_fasta (nifs, './tu_arch_seqs/tu_arch_')
