
###############################################################################
# LOADING LIBRARIES
###############################################################################


# How the library is structured (at present!)
import gene_cluster_library as gcl
import gene_cluster_analysis as gca
import gene_cluster_visualization as gcv


# Load the Stata nif library data
nifs = gcl.GeneClusterLibrary()
nifs.load('./data/nif_stata_library.txt')


# Can see basic structure with:
# nifs.parts
# nifs.variants


###############################################################################
# QUERYING LIBRARIES
###############################################################################


# Find instances of specific parts or types of part
p_insts = nifs.find_part_type_instances('Promoter')
p3_insts = nifs.find_part_instances('P3')


# Extract their sequence ranges (positions)
p_seq_ranges = nifs.find_seq_idx_ranges(p_insts)
p3_seq_ranges = nifs.find_seq_idx_ranges(p3_insts)


# Extract sequence with padding around part (useful for potential contextual effects)
p_seq_ranges_around = nifs.extract_seq_ranges(p_insts, 100, 50) # (instances, start offset, end offset)
p3_seq_ranges_around = nifs.extract_seq_ranges(p3_insts, 100, 50) # (instances, start offset, end offset)


# Extract all the divergent promoters
div_p = nifs.divergent_promoters()


# Extract all the divergent promoters
con_p = nifs.convergent_promoters()


# Extract monocistronic transcriptional units
monocis_tu = nifs.monocistronic_units()


# Extract polycistronic transcriptional units (with specific number of CDSs if needed)
polycis_tu = nifs.polycistronic_units(number_of_CDS=None)


# Scan for next types of parts (useful for query what is driven by promoters, etc)
p3_driven_cds = nifs.find_next_part_idxs(p3_insts, part_type='CDS', next_count=1, remove_nones=False, dir_check=False)


###############################################################################
# PLOTTING LIBRARIES
###############################################################################


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


if False:
	# Example plot testing full variant rendering
	fig = plt.figure(figsize=(14,3))
	ax = fig.add_subplot(3,1,1)
	gcv.plot_variant_arch(ax, nifs, '25', start_idx=1, end_idx=-1, linewidth=1.2)
	ax = fig.add_subplot(3,1,2)
	gcv.plot_variant_arch(ax, nifs, '10', start_idx=1, end_idx=-1, linewidth=1.2)
	ax = fig.add_subplot(3,1,3)
	gcv.plot_variant_arch(ax, nifs, '1', start_idx=1, end_idx=-1, linewidth=1.2)
	plt.tight_layout()


# Color mapping to use in library rendering
cmap = {}

# Promoters
cmap['P1'] = (1.0,0.0,0.0)
cmap['P2'] = (1.0,0.5,0.5)
cmap['P3'] = (1.0,0.8,0.8)

# Terminator
cmap['T1'] = (0.0,0.0,1.0)

# nif Genes
cmap['nifM'] = (0.5,0.77,0.56)
cmap['nifS'] = (0.37,0.78,0.78)
cmap['nifU'] = (0.76,0.47,0.76)
cmap['nifV'] = (0.78,0.78,0.49)
cmap['nifW'] = (0.94,0.74,0.4)
cmap['nifZ'] = (0.91,0.6,0.58)

# RBSs
rbsStrong = (0.16,0.68,0.15)
rbsWeak = (0.47,0.83,0.46)
cmap['Rm1'] = rbsStrong
cmap['Rm2'] = rbsWeak
cmap['Rs1'] = rbsStrong
cmap['Rs2'] = rbsWeak
cmap['Ru1'] = rbsStrong
cmap['Ru2'] = rbsWeak
cmap['Rv1'] = rbsStrong
cmap['Rv2'] = rbsWeak
cmap['Rw1'] = rbsStrong
cmap['Rw2'] = rbsWeak
cmap['Rz1'] = rbsStrong
cmap['Rz2'] = rbsWeak

# Hatch mapping to use in library rendering
hmap = {}
hmap['nifM'] = ''
hmap['nifS'] = '/'
hmap['nifU'] = '//'
hmap['nifV'] = '///'
hmap['nifW'] = '////'
hmap['nifZ'] = '/////'


if False:
	# Example plot of the whole library
	fig = plt.figure(figsize=(14,84))
	gcv.plot_library_arch(fig, nifs, linewidth=1.2, colormap=cmap, hatchmap=hmap)
	plt.tight_layout()


if True:
	# Example plot of architecture with Tx traces
	gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
	fig = plt.figure(figsize=(14,4))
	ax_arch = plt.subplot(gs[1])
	ax_traces = plt.subplot(gs[0],sharex=ax_arch)

	# Load the traces for predicted and measured
	phys_reads = gca.load_strand_data('./data/phys_depths3.csv')

	# Plot trace and the architecture
	gcv.plot_traces_with_arch(ax_arch, ax_traces, nifs, '79', phys_reads, start_idx=1, 
		                      end_idx=-1, linewidth=1.2, colormap=cmap, hatchmap=hmap)
	plt.tight_layout()


plt.show()

