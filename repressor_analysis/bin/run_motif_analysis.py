#!/usr/bin/env python
"""
Repressor analysis (E.coli DH10B)
"""

import motif_find as mf

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT\n\
			   Jing Zhang <jgzhang@mit.edu>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

###############################################################################
# Parameters
genome_filename = '../data/genomes/DH10B.fasta'
repressors = ['AmeR', 'AmtR', 'ButR', 'LitR', 'McbR', 'Orf2', 
              'PsrA', 'QacR', 'ScbR', 'SrpR', 'TarA']
###############################################################################

def run_analysis (repressor_filename, genome_seq, background, output_prefix):
	"""All analysis is coordinated from here
	"""
	# Load the motif
	print('\tLoading motifs')
	motif = mf.load_motifs(repressor_filename)
	# Normalise to the approx base distribution of the host (avoid overfitting)
	print('\tGenerating PWM and PSSM')
	pwm = motif.counts.normalize(pseudocounts=0.5)
	pssm = pwm.log_odds(background)
	# Calculate threshold to use
	print('\tCalculating threshold')
	distribution = pssm.distribution(background=background, precision=10**4)
	threshold = distribution.threshold_patser()
	print('\tThreshold = '+str(threshold))
	# Find hits
	print('\tSearching for score hits above threshold')
	score_hits = mf.search_score_hits(pssm, genome_seq, threshold=threshold)
	# Save data to file
	print('\tSaving score hits')
	f_out = open(output_prefix+'score_hits.csv', 'w')
	for c in score_hits.keys():
		f_out.write('>'+c+'\n')
		for el in score_hits[c]:
			f_out.write(','.join([str(x) for x in el])+'\n')
	f_out.close()

# Load data files
genome_seq = mf.load_genome_seq(genome_filename)
background = mf.genome_background(genome_seq)

# Process each repressor
for r in repressors:
	print('Processing '+r+'...')
	run_analysis('../data/motifs/'+r+'.txt', genome_seq, 
		         background, '../results/'+r+'_')
print('Done.')
