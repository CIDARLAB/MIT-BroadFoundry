#!/usr/bin/env python
"""
Motif Find
===========
	Will search for potential binding sites of repressors based on a PWM matrix
	of the binding seqeunces. See main method for details on calling function.
"""

import csv
import numpy as np
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

def load_metadata (metadata_filename):
	metadata = {}
	with open(metadata_filename, 'rU') as metadatafile:
		metadatareader = csv.reader(metadatafile, delimiter='\t')
		# Ignore header
		metadatareader.next()
		for row in metadatareader:
			metadata[row[0]] = row[1:]
	return metadata

def load_tss (tss_filename):
	tss = {}
	with open(tss_filename, 'rU') as tssfile:
		tssreader = csv.reader(tssfile, delimiter='\t')
		# Ignore header
		tssreader.next()
		for row in tssreader:
			if row[0] not in tss.keys():
				tss[row[0]] = {}

			# TSS positions are in 1...len (convert to index)
			tss[row[0]][row[1]] = [[x.strip() for x in row[2].split(',')],
			                       row[3],
			                       int(row[4])-1,
			                       int(row[5])]
	return tss

def load_fpkms (fpkm_filename):
	fpkms = {}
	with open(fpkm_filename, 'rU') as fpkmfile:
		fpkmreader = csv.reader(fpkmfile, delimiter='\t')
		# Ignore header
		fpkmreader.next()
		for row in fpkmreader:
			if row != '':
				if row[0] not in fpkms.keys():
					fpkms[row[0]] = {}
				fpkms[row[0]][row[1]] = [float(x) for x in row[3:]]
	return fpkms

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

def search_pssm_hits (pssm, genome_seq, threshold=0.0, relative=False):
	"""Search for positions where score is above given threshold
	   Positions are given as +ve and -ve indexes to use in: 
	   seq[pos:pos+len(motif)]
	"""
	if relative == True:
		threshold = (pssm.max - pssm.min) * threshold + pssm.min
	scores = {}
	for c in genome_seq.keys():
		scores[c] = []
		for position, score in pssm.search(genome_seq[c], threshold=threshold):
			rel_score = (score - pssm.min) / (pssm.max - pssm.min)
			scores[c].append([position, score, rel_score])
	return scores

def total_hits (scores):
	count = 0
	for c in scores.keys():
		count += len(scores[c])
	return count

def load_score_hits (score_filename):
	scores = {}
	cur_chr = ''
	with open(score_filename, 'rU') as scorefile:
		scorereader = csv.reader(scorefile, delimiter=',')
		for row in scorereader:
			if row[0] != '':
				if row[0][0] == '>':
					cur_chr = row[0][1:]
					if cur_chr in scores.keys():
						scores[cur_chr] = np.array(scores[cur_chr])
					scores[cur_chr] = []
				else:
					scores[cur_chr].append([int(row[0]), float(row[1])])
	scores[cur_chr] = np.array(scores[cur_chr])
	return scores

def extract_tss_hits (motif, genome_len, tss, scores, us_window, ds_window):
	tss_hits = {}
	for c in tss.keys():
		for tss_key in tss[c].keys():
			cur_tss = tss[c][tss_key]
			cur_genes = cur_tss[0]
			cur_dir = cur_tss[1]
			cur_pos = cur_tss[2]
			cur_atg_pos = [3]
			start_pos = 0
			end_pos = 0
			if c not in tss_hits.keys():
				tss_hits[c] = []
			if cur_dir == '+':
				start_pos = cur_pos-us_window
				end_pos = cur_pos+ds_window
			else:
				start_pos = cur_pos+us_window
				end_pos = cur_pos-ds_window
			for hit in scores[c]:
				hit_pos = hit[0]
				if hit_pos < 0:
					hit_pos = genome_len+hit_pos
				if hit_pos+len(motif) >= start_pos and hit_pos <= end_pos:
					tss_hits[c].append(hit)
	for c in tss_hits.keys():
		tss_hits[c] = np.array(tss_hits[c])
	return tss_hits

def filter_hits (tss_hits, threshold):
	filtered_tss_hit = {}
	for c in tss_hits.keys():
		for hit_idx in range(np.size(tss_hits[c], axis=0)):
			cur_hit = list(tss_hits[c][hit_idx,:])
			if cur_hit[1] >= threshold:
				if c not in filtered_tss_hit.keys():
					filtered_tss_hit[c] = []
				filtered_tss_hit[c].append(cur_hit)
	for c in filtered_tss_hit.keys():
		filtered_tss_hit[c] = np.array(filtered_tss_hit[c])
	return filtered_tss_hit

def hit_per_tss (tss, tss_hits, motif, genome_len, us_window, ds_window):
	hits_per_tss = {}
	for c in tss.keys():
		hits_per_tss[c] = {}
		for tss_key in tss[c].keys():
			cur_tss = tss[c][tss_key]
			cur_genes = cur_tss[0]
			cur_dir = cur_tss[1]
			cur_pos = cur_tss[2]
			cur_atg_pos = [3]
			start_pos = 0
			end_pos = 0
			if cur_dir == '+':
				start_pos = cur_pos-us_window
				end_pos = cur_pos+ds_window
			else:
				start_pos = cur_pos+us_window
				end_pos = cur_pos-ds_window

			# Go through and extract hits above threshold
			for hit in tss_hits[c]:
				hit_pos = hit[0]
				if hit_pos < 0:
					hit_pos = genome_len+hit_pos
				
				if hit_pos-len(motif) >= start_pos and hit_pos <= end_pos:
					if cur_dir == '+':
						bind_start = (hit_pos-cur_pos)
						bind_end = (hit_pos-cur_pos)+len(motif)
					else:
						bind_start = (cur_pos-hit_pos)
						bind_end = (cur_pos-hit_pos)+len(motif)
					if tss_key not in hits_per_tss[c].keys():
						hits_per_tss[c][tss_key] = []
					hits_per_tss[c][tss_key].append([bind_start, bind_end, hit[1]])

	return hits_per_tss

def binding_trace (tss, tss_hits, motif, genome_len, us_window, ds_window, threshold):
	x_vals = np.array(range(-us_window,ds_window+1))
	bindings = np.zeros(len(x_vals))

	for c in tss.keys():
		for tss_key in tss[c].keys():
			cur_tss = tss[c][tss_key]
			cur_genes = cur_tss[0]
			cur_dir = cur_tss[1]
			cur_pos = cur_tss[2]
			cur_atg_pos = [3]
			start_pos = 0
			end_pos = 0
			if cur_dir == '+':
				start_pos = cur_pos-us_window
				end_pos = cur_pos+ds_window
			else:
				start_pos = cur_pos+us_window
				end_pos = cur_pos-ds_window

			# Go through and extract hits above threshold
			for hit in tss_hits[c]:
				hit_pos = hit[0]
				if hit_pos < 0:
					hit_pos = genome_len+hit_pos
				
				if hit[1] >= threshold and hit_pos+len(motif) >= start_pos and hit_pos <= end_pos:
					# Find hit position and increment bindings over motif range
					bind_start = 0
					bind_end = 0

					if cur_dir == '+':
						bind_start = us_window+(hit_pos-cur_pos)
						if bind_start < 0:
							bind_start = 0
						bind_end = us_window+(hit_pos-cur_pos)+len(motif)
						if bind_end > len(x_vals):
							bind_end = len(x_vals)
					else:
						bind_start = us_window+(cur_pos-hit_pos)
						if bind_start < 0:
							bind_start = 0
						bind_end = us_window+(cur_pos-hit_pos)+len(motif)
						if bind_end > len(x_vals):
							bind_end = len(x_vals)

					for idx in range(int(bind_start), int(bind_end)):
						bindings[idx] += 1

	return x_vals, bindings
