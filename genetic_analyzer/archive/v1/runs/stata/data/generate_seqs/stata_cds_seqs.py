#!/usr/bin/env python
"""Extract all CDSs from a multi-FASTA and GFF annotation
"""

import csv
from Bio import SeqIO

def reverse_complement(seq):
	"""http://crazyhottommy.blogspot.com/2013/10/python-code-for-getting-reverse.html
	"""
	for base in seq:
		if base not in 'ATCGatcg':
			print "Error: NOT a DNA sequence"
			return None
	seq1 = 'ATCGTAGCatcgtagc'
	seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
	return "".join([seq_dict[base] for base in reversed(seq)])

# Load the multi-FASTA file to dict on contig name
contigs = {}
contig_seqs = SeqIO.parse(open('./refs/Klebsiella_oxytoca_M5al.fasta'),'fasta')
for contig in contig_seqs:
	name = contig.id
	sequence = contig.seq.tostring()
	contigs[name.split('|')[3]] = sequence

# Parse GFF file and extract seq for each gene (save to multi-FASTA output)
gene_seqs = {}
gff_reader = csv.reader(open('./refs/Klebsiella_oxytoca_M5al.gff', 'rU'), delimiter='\t')
for row in gff_reader:
	if row >= 9:
		gene_name = row[8].split('=')[1]
		start_idx = int(row[3])-1
		end_idx = int(row[4])
		strand = row[6]
		contig = row[0]
		# Extract the gene sequence
		contig_seq = contigs[contig]
		gene_seq = contig_seq[start_idx:end_idx].upper()
		# Check strand and reverse complement if necessary
		if strand == '-':
			gene_seq = reverse_complement(gene_seq)
		gene_seqs[gene_name] = gene_seq

# Generate gene sequences for host
fasta_out = open('../host_gene_seqs.fa', 'w')
for gene in sorted(gene_seqs.keys()):
	fasta_out.write('>'+gene+'\n'+gene_seqs[gene]+'\n\n')
fasta_out.close()

# Load the multi-FASTA file to dict on contig name
contigs = {}
contig_seqs = SeqIO.parse(open('./refs/synnifI4.fasta'),'fasta')
for contig in contig_seqs:
	name = contig.id
	sequence = contig.seq.tostring().upper()
	contigs[name] = sequence

# Parse GFF file and extract seq for each gene (save to multi-FASTA output)
gene_seqs = {}
gff_reader = csv.reader(open('./refs/synnifI4.gff', 'rU'), delimiter='\t')
for row in gff_reader:
	if row >= 9:
		gene_name = row[8].split('=')[1]
		start_idx = int(row[3])-1
		end_idx = int(row[4])
		strand = row[6]
		contig = row[0]
		# Extract the gene sequence
		contig_seq = contigs[contig]
		gene_seq = contig_seq[start_idx:end_idx].upper()
		# Check strand and reverse complement if necessary
		if strand == '-':
			gene_seq = reverse_complement(gene_seq)
		gene_seqs[gene_name] = gene_seq

# Generate gene sequences for synthetic genes
fasta_out = open('../syn_gene_seqs.fa', 'w')
for gene in sorted(gene_seqs.keys()):
	fasta_out.write('>'+gene+'\n'+gene_seqs[gene]+'\n\n')
fasta_out.close()
