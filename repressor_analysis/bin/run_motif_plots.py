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
tss_filename = '../data/tss/DH10B_TSS.tsv'
us_window = 200
ds_window = 200
###############################################################################

def run_analysis (repressor_filename, genome_seq, background, tss, output_prefix):
	"""All analysis is coordinated from here
	"""
	# Load the motif
	motif = mf.load_motifs(repressor_filename)
	# Normalise to the approx base distribution of the host (avoid overfitting)
	pwm = motif.counts.normalize(pseudocounts=0.5)
	pssm = pwm.log_odds(background)
	# Generate score profile
	full_scores = mf.calculate_scores(pssm, genome_seq)
	# Calculate threshold to use
	#distribution = pssm.distribution(background=background, precision=10**4)
	#threshold = distribution.threshold_patser()
	threshold = 0.0
	# Find hits
	scores = mf.search_score_hits(pssm, genome_seq, threshold=threshold)
	hit_count = mf.total_hits(scores)
	tss_hits = mf.extract_tss_hits (tss, scores, us_window, ds_window)

# Load data files
genome_seq = mf.load_genome_seq(genome_filename)
background = mf.genome_background(genome_seq)
#tss = mf.load_tss(tss_filename)
tss = None

# Process each repressor
for r in ['AmeR', 'AmtR', 'ButR', 'LitR', 'McbR', 'Orf2', 
          'PsrA', 'QacR', 'ScbR', 'SrpR', 'TarA']:
	print("Processing "+r)
	run_analysis('../data/motifs/'+r+'.txt', genome_seq, 
		         background, tss, '../results/'+r+'_')
