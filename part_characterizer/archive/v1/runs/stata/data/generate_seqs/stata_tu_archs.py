
import subprocess
import gene_cluster_library as gct
import gene_cluster_analysis as gca

# Load the example data set
nifs = gct.GeneClusterLibrary()
nifs.load('../stata_library.txt')

# Extract all the sequences for the TUs (we can use this data for fitting)
gca.generate_all_variant_tus_fasta(nifs, '../tu_archs/tu_arch_', non_terminated=True, 
	                               truncate_ends=1, none_end_idx_offset=0, 
	                               none_end_bp_pad=80, none_start_bp_pad=80, 
	                               id_prefix='SYNTHETIC_')
