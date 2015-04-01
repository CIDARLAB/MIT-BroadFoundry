#!/usr/bin/env python
"""
Motif Find
===========
	Will search for potential binding sites of repressors based on a PWM matrix
	of the binding seqeunces. See main method for details on calling function.
"""

import csv
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA

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
		name, sequence = fasta.id, str(fasta.seq)
		genome_seqs[name] = Seq(sequence, IUPACUnambiguousDNA())
	return genome_seqs

def load_tss (tss_filename):
	tss = {}
	with open('tss_filename.csv', 'rU') as tssfile:
		tssreader = csv.reader(tssfile, delimiter='\t')
		for row in tssreader:
			tss[row[0]] = [row[1], 
			               [x.strip() for x in row[2].split(',')],
			               row[3],
			               int(row[4])]
	return tss

def genome_background (genome_seq):
	"""Generate a background fraction for each base from genome
	"""
	background = {}
	counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
	total_count = 0.0
	for c in genome_seq.keys():
		for b in 'ATGC':
			counts[b] += genome_seq[c].count(b)
		total_count += len(genome_seq[c])
	for b in 'ATGC':
		background[b] = float(counts[b])/total_count
	return background

def calculate_scores (pssm, genome_seq):
	"""Calculate the scores for a matrix on +ve and -ve strands
	"""
	scores = {}
	rpssm = pssm.reverse_complement()
	for c in genome_seq.keys():
		scores[c] = {}
		scores[c]['+'] = pssm.calculate(genome_seq[c])
		scores[c]['-'] = rpssm.calculate(genome_seq[c])
	return scores

def search_score_hits (pssm, genome_seq, threshold=0.0):
	"""Search for positions where score is above given threshold
	   Positions are given as +ve and -ve indexes to use in: 
	   seq[pos:pos+len(motif)]
	"""
	scores = {}
	for c in genome_seq.keys():
		scores[c] = []
		for position, score in pssm.search(genome_seq[c], threshold=threshold):
			scores[c].append([position, score])
	return scores

def total_hits (scores):
	count = 0
	for c in scores.keys():
		count += len(scores[c])
	return count

def extract_tss_hits (tss, scores, us_window, ds_window):
	
	return None
