"""
	Plot part strengths from RSEM generated promoter and terminator performance values
"""

# Required modules
import csv
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

###############################################################################
# PARAMETERS
###############################################################################

# Which data to process?
# STATA LIBRARY
#LIBRARY_TO_PROCESS = 'stata_trunc'
# CIRCUIT WORKING/BROKEN
LIBRARY_TO_PROCESS = 'circuit'

# Measure to use to calculate performance
MEASURE_USED = 'fpkm'

# Standard locations
RESULTS_PREFIX = '../results/'+LIBRARY_TO_PROCESS+'/part_perf/'

###############################################################################
# SUPPORTING FUNCTIONS
###############################################################################

def load_perf_data (filename):
	perf_data = {}
	file_reader = csv.reader(open(filename, 'rU'), delimiter=',')
	# Ignore the header
	header = next(file_reader)
	# Process each line
	for row in file_reader:
		if len(row) > 3:
			variant = row[0]
			part_idx = int(row[1])
			part_name = row[2]
			tx_data = [float(x) for x in row[3:-1]]
			strength = float(row[-1])
			if part_name not in perf_data.keys():
				perf_data[part_name] = []
			perf_data[part_name].append([variant, part_idx, tx_data, strength])
	return perf_data

###############################################################################
# PLOT THE RESULTS
###############################################################################

# Load the data
p_perf_data = load_perf_data(RESULTS_PREFIX+'promoter_express_pref.txt')
t_perf_data = load_perf_data(RESULTS_PREFIX+'terminator_express_pref.txt')
# Plot the promoter distribution
fig = plt.figure(figsize=(9,6))
x_min = 0
x_max = 1000000
bin_size = 50000
p_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
ax_idx = 1
for promoter in sorted(p_perf_data.keys()):
	# Extract the data
	p_data = [x[-1] for x in p_perf_data[promoter]]
	# Plot the data
	ax = plt.subplot(2,len(p_perf_data.keys())/2,ax_idx)
	ax.hist(p_data, bins=p_bin_list, color='blue', alpha=0.2)
	plt.axvline(x=np.median(p_data), color='blue', linestyle='-', linewidth=1.5)
	ax.set_title(promoter + '('+str(np.median(p_data))+')')
	ax.set_xlabel('Strength (FPKM)')
	ax.set_ylabel('Count')
	ax_idx += 1
plt.tight_layout()
plt.savefig('./rsem_stata_promoter_perf.pdf', transparent=True)
plt.close('all')

# Plot the promoter in -> out fluxes
fig = plt.figure(figsize=(9,6))
ax_idx = 1
for promoter in sorted(p_perf_data.keys()):
	# Extract the data
	p_tx_in = [x[2][0] for x in p_perf_data[promoter]]
	p_tx_out = [x[-1] for x in p_perf_data[promoter]]
	# Plot the data
	ax = plt.subplot(2,len(p_perf_data.keys())/2,ax_idx)
	ax.scatter(p_tx_in, p_tx_out)
	ax.set_title(promoter)
	ax.set_xlabel('Tx IN (FPKM)')
	ax.set_ylabel('Strength (FPKM)')
	ax_idx += 1
plt.tight_layout()
plt.savefig('./rsem_stata_promoter_tx_in_out.pdf', transparent=True)
plt.close('all')

# Plot the terminator distribution
fig = plt.figure(figsize=(9,6))
x_min = 0
x_max = 1
bin_size = 0.05
t_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
ax_idx = 1
for terminator in sorted(t_perf_data.keys()):
	# Extract the data
	t_data = [x[-1] for x in t_perf_data[terminator]]
	# Plot the data
	ax = plt.subplot(2,len(t_perf_data.keys())/2,ax_idx)
	ax.hist(t_data, bins=t_bin_list, color='blue', alpha=0.2)
	plt.axvline(x=np.median(t_data), color='blue', linestyle='-', linewidth=1.5)
	ax.set_title(terminator + '('+str(np.median(t_data))+')')
	ax.set_xlabel('Strength (% FPKM)')
	ax.set_ylabel('Count')
	ax_idx += 1
plt.tight_layout()
plt.savefig('./rsem_stata_terminator_perf.pdf', transparent=True)
plt.close('all')

# Plot the terminator in -> out fluxes
fig = plt.figure(figsize=(9,6))
x_min = 0
x_max = 1
bin_size = 0.05
t_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
ax_idx = 1
for terminator in sorted(t_perf_data.keys()):
	# Extract the data
	t_tx_in = [x[2][0]/10000.0 for x in t_perf_data[terminator] if x[0][2] == '2' and x[0][0] in ['6','7','8']] #and x[0][0] in ['6','7','8'] 
	t_tx_out = [x[-1] for x in t_perf_data[terminator] if x[0][2] == '2' and x[0][0] in ['6','7','8']]
	#print t_tx_in
	#print t_tx_in
	# Plot the data
	ax = plt.subplot(2,len(t_perf_data.keys())/2,ax_idx)
	ax.scatter(t_tx_in, t_tx_out)
	ax.set_title(terminator)
	ax.set_xlabel('In Flux (FPKM)')
	ax.set_ylabel('% Termination (FPKM)')
	ax.set_ylim([0,1.05])
	#ax.set_xlim([0.01, np.max(t_tx_in)*1.05])
	ax.set_xlim([0.1 , 20])

	plt.axhline(y=np.median(t_tx_out), color='red')
	ax.text(5, 0.55, str(np.median(t_tx_out)*100)[0:4])

	#xi = np.arange(0.0, np.max(t_tx_in)*1.05, )
	#slope, intercept, r_value, p_value, std_err = stats.linregress(t_tx_in,t_tx_out)
	#line = slope*xi+intercept
	#ax.plot(xi,line,'r-')
	#ax.text(0, 0.8, 'R: ' + str(r_value))
	#ax.text(0, 0.7, 'a: ' + str(slope))
	#ax.text(0, 0.6, 'p: ' + str(p_value))
	ax.set_xscale('log')
	#ax.set_yscale('log')
	ax_idx += 1
plt.tight_layout()
plt.savefig('./WFAILURE_rsem_stata_terminator_tx_in_out.pdf', transparent=True)
plt.close('all')




# Plot the terminator in -> out fluxes
fig = plt.figure(figsize=(9,6))
x_min = 0
x_max = 1
bin_size = 0.05
t_bin_list = np.arange(x_min, x_max + bin_size, bin_size)
ax_idx = 1
for terminator in sorted(t_perf_data.keys()):
	# Extract the data
	t_tx_in = [x[2][0]/10000.0 for x in t_perf_data[terminator] if x[0][2] == '1' and x[0][0] in ['6','7','8']] #and x[0][0] in ['6','7','8'] 
	t_tx_out = [x[-1] for x in t_perf_data[terminator] if x[0][2] == '1' and x[0][0] in ['6','7','8']]
	#print t_tx_in
	#print t_tx_in
	# Plot the data
	ax = plt.subplot(2,len(t_perf_data.keys())/2,ax_idx)
	ax.scatter(t_tx_in, t_tx_out)
	ax.set_title(terminator)
	ax.set_xlabel('In Flux (FPKM)')
	ax.set_ylabel('% Termination (FPKM)')
	ax.set_ylim([0,1.05])
	#ax.set_xlim([0.01, np.max(t_tx_in)*1.05])
	ax.set_xlim([0.1 , 20])

	plt.axhline(y=np.median(t_tx_out), color='red')
	ax.text(5, 0.55, str(np.median(t_tx_out)*100)[0:4])

	#xi = np.arange(0.0, np.max(t_tx_in)*1.05, )
	#slope, intercept, r_value, p_value, std_err = stats.linregress(t_tx_in,t_tx_out)
	#line = slope*xi+intercept
	#ax.plot(xi,line,'r-')
	#ax.text(0, 0.8, 'R: ' + str(r_value))
	#ax.text(0, 0.7, 'a: ' + str(slope))
	#ax.text(0, 0.6, 'p: ' + str(p_value))
	ax.set_xscale('log')
	#ax.set_yscale('log')
	ax_idx += 1
plt.tight_layout()
plt.savefig('./WORKING_rsem_stata_terminator_tx_in_out.pdf', transparent=True)
plt.close('all')