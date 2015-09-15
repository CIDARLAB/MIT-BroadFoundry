#!/usr/bin/env python
"""Extract all CDSs from a Genbank files for host and circuit
"""

from Bio import SeqIO

# Generate gene sequences for synthetic circuit
fasta_out = open('../syn_gene_seqs.fa', 'w')
for rec in SeqIO.parse('pAAKN-0x58v50.gb', 'genbank'):
	if rec.features:
		for feature in rec.features:
			if feature.type == 'misc_feature':
				cds_name = feature.qualifiers['label'][0] + \
				           '|' + feature.qualifiers['label'][0] + \
				           '|' + str(feature.location)
				cds_seq = feature.location.extract(rec).seq
				# Write the record to file
				fasta_out.write('>SYNTHETIC_'+cds_name+'\n'+str(cds_seq)+'\n\n')
fasta_out.close()

extracted_features = 0

# Generate gene sequences for host
fasta_out = open('../host_gene_seqs.fa', 'w')
for rec in SeqIO.parse('Escherichia_coli_DH10B.gb', 'genbank'):
	if rec.features:
		for feature in rec.features:
			if feature.type == 'gene': #CDS
				cds_name = feature.qualifiers['locus_tag'][0]
				if 'gene' in feature.qualifiers.keys():
					cds_name = cds_name + '|' + feature.qualifiers['gene'][0]
				else:
					cds_name = cds_name + '|' + 'no_gene_name'
				cds_name = cds_name + '|' + str(feature.location)
				cds_seq = feature.location.extract(rec).seq
				# Write the record to file
				fasta_out.write('>'+cds_name+'\n'+str(cds_seq)+'\n\n')
				extracted_features +=1
fasta_out.close()

print 'Extracted:', extracted_features, 'features'
