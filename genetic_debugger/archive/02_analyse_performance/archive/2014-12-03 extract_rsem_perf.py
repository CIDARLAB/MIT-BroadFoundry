"""
	From an RSEM run extract the promoter and terminator performance values
"""

# Required modules
import csv
import numpy as np
import gene_cluster_library as gcl

###############################################################################
# PARAMETERS
###############################################################################

LIBRARY_TO_PROCESS = 'stata_trunc'
#LIBRARY_TO_PROCESS = 'stata'

GCL_DATA = '../data/'+LIBRARY_TO_PROCESS+'/clean_nif_stata_library.txt'
RSEM_DATA_PREFIX = '../results/'+LIBRARY_TO_PROCESS+'-ss-pe/rsem_exps/'
MEASURE_TO_USE = 'TPM'
RESULTS_PREFIX = '../results/'+LIBRARY_TO_PROCESS+'-ss-pe/rsem_perf/'

# Which designs to process
designs_to_process = [str(x) for x in range(1,85)]
# Remove designs with sequencing errors
temp_designs = []
for idx in range(len(designs_to_process)):
	if designs_to_process[idx] not in ['17', '18', '33', '45', '57', '69', '76', '81']:
		temp_designs.append(designs_to_process[idx])
designs_to_process = temp_designs

###############################################################################
# SUPPORTING FUNCTIONS
###############################################################################

def load_rsem_gene_results (filename):
	"""Load the RSEM gene results file into a dict of dicts.
	"""
	rsem_data = []
	file_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Process header for column names (used in dict returned)
	header = next(file_reader)
	header_map = {}
	for idx in range(len(header)):
		header_map[idx] = header[idx]
	# Process each line
	for row in file_reader:
		if len(row) == 7:
			tx_data = {}
			tx_data[header_map[2]] = float(row[2])
			tx_data[header_map[3]] = float(row[3])
			tx_data[header_map[4]] = float(row[4])
			tx_data[header_map[5]] = float(row[5])
			tx_data[header_map[6]] = float(row[6])
			# Extract other attributes
			id_parts = row[0].split('_')
			tx_data['dir'] = id_parts[1]
			p_idxs = id_parts[2].split('-')
			p_bps = id_parts[3].split('-')
			if p_idxs[0] == 'None':
				tx_data['start_idx'] = None
				tx_data['start_bp'] = None
			else:
				tx_data['start_idx'] = int(p_idxs[0])
				tx_data['start_bp'] = int(p_bps[0])
			if p_idxs[1] == 'None':
				tx_data['end_idx'] = None
				tx_data['end_bp'] = None
			else:
				tx_data['end_idx'] = int(p_idxs[1])
				tx_data['end_bp'] = int(p_bps[1])
		# Add the data element
		rsem_data.append(tx_data)
	return rsem_data

###############################################################################
# EXTRACT THE ESTIMATES
###############################################################################

# Load the design data (needed to link promoter/terminator names to part idxs)
gcl_data = gcl.GeneClusterLibrary()
gcl_data.load(GCL_DATA)

# -----------------
# PROMOTER ANALYSIS
# -----------------

# List of the extracted performances
p_perf_data = {}

# For each design extract performance data
for design in designs_to_process:
	# Add placeholder for data
	p_perf_data[design] = {}

	# Load the data
	design_file = RSEM_DATA_PREFIX+'rsem_exp_'+design+'.genes.results'
	design_data = load_rsem_gene_results(design_file)

	# For each transript 
	for trans in design_data:
		# Find promoter 
		promoter_idx = trans['start_idx']

		# Have we seen this promoter idx before?
		if promoter_idx not in p_perf_data[design].keys():
			p_perf_data[design][promoter_idx] = []

		# Add the performance data
		p_perf_data[design][promoter_idx].append(trans[MEASURE_TO_USE])

# Write the performance data to file
p_out = open(RESULTS_PREFIX+'promoter_rsem_pref.txt', 'w')
# Write header
p_out.write('variant,promoter_idx,promoter_name,strength('+MEASURE_TO_USE+')\n')

TEMP_P1 = []
TEMP_P2 = []
TEMP_P3 = []
for design in sorted(p_perf_data.keys()):
	for promoter_idx in p_perf_data[design]:
		p_strength = np.sum(p_perf_data[design][promoter_idx])
		data_out = [design, str(promoter_idx), gcl_data.variant_part_idx_name(design, promoter_idx),
		            str(p_strength)+'\n']
		p_out.write(','.join(data_out))

		if gcl_data.variant_part_idx_name(design, promoter_idx) == 'P1':
			TEMP_P1.append(p_strength)
		if gcl_data.variant_part_idx_name(design, promoter_idx) == 'P2':
			TEMP_P2.append(p_strength)
		if gcl_data.variant_part_idx_name(design, promoter_idx) == 'P3':
			TEMP_P3.append(p_strength)
p_out.close()

# -------------------
# TERMINATOR ANALYSIS
# -------------------

# List of the extracted performances
t_perf_data = {}

# For each design extract performance data
for design in designs_to_process:
	# Add placeholder for data
	t_perf_data[design] = {}

	# Load the data
	design_file = RSEM_DATA_PREFIX+'rsem_exp_'+design+'.genes.results'
	design_data = load_rsem_gene_results(design_file)

	# Find all the terminators in the transcripts (check type of end_idx)
	t_list = []
	for trans in design_data:
		part_idx = trans['end_idx']
		if part_idx != None and gcl_data.variant_part_idx_type(design, part_idx) == 'Terminator':
				if part_idx not in t_list:
					t_list.append(part_idx)

	# Process each terminator
	for term_idx in t_list:
		term_dir = gcl_data.variant_part_idx_dir(design, term_idx)
		tx_before = []
		tx_after = []
		# Go through each transcript and classify if in before or after group
		for trans in design_data:
			if term_dir == 'F':
				if trans['start_idx'] < term_idx and trans['end_idx'] == term_idx:
					tx_before.append(trans[MEASURE_TO_USE])
				if trans['start_idx'] < term_idx and trans['end_idx'] > term_idx:
					tx_after.append(trans[MEASURE_TO_USE])
			else:
				if trans['start_idx'] > term_idx and trans['end_idx'] == term_idx:
					tx_before.append(trans[MEASURE_TO_USE])
				if trans['start_idx'] > term_idx and trans['end_idx'] < term_idx:
					tx_after.append(trans[MEASURE_TO_USE])
		if np.sum(tx_before) == 0.0 and np.sum(tx_after) == 0.0:
			t_perf_data[design][term_idx] = 0.0
		else:
			t_perf_data[design][term_idx] = np.sum(tx_before)/(np.sum(tx_before)+np.sum(tx_after))

# Write the performance data to file
t_out = open(RESULTS_PREFIX+'terminator_rsem_pref.txt', 'w')
# Write header
t_out.write('variant,terminator_idx,terminator_name,strength('+MEASURE_TO_USE+')\n')

TEMP_T = []
for design in sorted(t_perf_data.keys()):
	for terminator_idx in t_perf_data[design]:
		data_out = [design, str(terminator_idx), gcl_data.variant_part_idx_name(design, terminator_idx),
		            str(t_perf_data[design][terminator_idx])+'\n']
		t_out.write(','.join(data_out))
		TEMP_T.append(t_perf_data[design][terminator_idx])
t_out.close()

# --------------------------
# PLOT OVERVIEW OF STRENGTHS
# --------------------------

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(5,5))

ax1 = plt.subplot(2,1,1)
x_min = 0
x_max = 1000000
bin_size = 50000
p_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
ax1.hist(TEMP_P1, bins=p_bin_list, color='blue', alpha=0.2, label='P1')
plt.axvline(x=np.median(TEMP_P1), color='blue', linestyle='-', linewidth=1.5)
ax1.hist(TEMP_P2, bins=p_bin_list, color='green', alpha=0.2, label='P2')
plt.axvline(x=np.median(TEMP_P2), color='green', linestyle='-', linewidth=1.5)
ax1.hist(TEMP_P3, bins=p_bin_list, color='red', alpha=0.2, label='P3')
plt.axvline(x=np.median(TEMP_P3), color='red', linestyle='-', linewidth=1.5)
ax1.legend(fontsize=12)
ax1.set_title('Promoter Strengths')
ax1.set_xlabel('Strength (TPM)')
ax1.set_ylabel('Count')


ax2 = plt.subplot(2,1,2)
x_min = 0
x_max = 1
bin_size = 0.05
term_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
ax2.hist(TEMP_T, bins=term_bin_list, color=[0.5,0.5,0.5], edgecolor=[0.5,0.5,0.5])
plt.axvline(x=np.median(TEMP_T), color=(0.2,0.2,0.2), linestyle='-', linewidth=1.5)
ax2.set_title('Terminator Efficiencies')
ax2.set_xlabel('Termination Efficiency (%)')
ax2.set_ylabel('Count')

plt.tight_layout()
plt.savefig('./rsem_perf_overview.pdf', transparent=True)
plt.close('all')
