#!/usr/bin/env python
"""
Motif Find
===========
	Will search for potential binding sites of repressors based on a PWM matrix
	of the binding seqeunces. See main method for details on calling function.
"""
#    Motif Find
#    Copyright (C) 2015 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    Jing Zhang <jgzhang@mit.edu>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

from argparse import ArgumentParser
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT\n\
			   Jing Zhang <jgzhang@mit.edu>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

def load_motifs (motif_filename):
	"""Load seqence motifs
	"""
	instances = []
	f = open(motif_filename, 'rU')
	for line in f:
		if line.strip() != '':
			instances.append(Seq(line.strip()))
	return motifs.create(instances)

def load_genome_seq (genome_filename):
	"""Load genome sequence from multi-FASTA (dict keyed on chromosome name)
	"""
	genome_seqs = {}
	fasta_sequences = SeqIO.parse(open(genome_filename),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        genome_seqs[name] = Seq(sequence)
	return genome_seqs



print(load_motifs('./data/AmeR.csv'))



instances = [Seq("TACAA"), Seq("TATTA"), Seq("GGGAA"), Seq("TACAC")]
m = motifs.create(instances)
print(m.counts)
print(m.consensus)
r = m.reverse_complement()
print(r.consensus)

# Normalise to the approx base distribution of the host (avoid overfitting)
pwm = m.counts.normalize(pseudocounts=0.5)
pwm = m.counts.normalize(pseudocounts={'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6})
print(pwm)

# Position-specific scoring matrix
background = {'A':0.3,'C':0.2,'G':0.2,'T':0.3}
pssm = pwm.log_odds(background)


test_seq=Seq("TACACTGCATTACAACCCAAGCATTA", m.alphabet)
len(test_seq)

# threshold of 0 is more likely to be motif than background.
for position, score in pssm.search(test_seq, threshold=3.0):
	print("Position %d: score = %5.3f" % (position, score))


# positions are are 
test_seq[pos:pos+len(m)]

# for all positions
pssm.calculate(test_seq)
rpssm = pssm.reverse_complement()
rpssm.calculate(test_seq)



# selectiong threshold
distribution = pssm.distribution(background=background, precision=10**4)

# based on FDR
threshold = distribution.threshold_fpr(0.01)













