#!/usr/bin/env python

#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Requires following software is available and in path:
#   - BWA
#   - SAMTools
#   - HTSeq
#   - BEDTools
#   - R + edgeR package (Rscript)

# Required modules
import csv
import subprocess

def bwa_index_filename (settings, sample):
	return settings[sample]['temp_path']+sample

def sam_filename (settings,  sample):
	return settings[sample]['temp_path']+sample+'.sam'

def bam_filename (settings,  sample, extension=True):
	bam_filename = settings[sample]['temp_path']+sample
	if extension == True:
		bam_filename += '.bam'
	return bam_filename

def count_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.counts.txt'

def mapped_reads_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.mapped.reads.txt'

def mapped_reads_matrix_filename (settings):
	return settings['None']['output_path']+'mapped.reads.matrix.txt'

def gene_length_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.gene.lengths.txt'

def profile_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.profiles.txt'

def profile_fwd_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fwd.profiles.txt'

def profile_rev_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.rev.profiles.txt'

def count_matrix_filename (settings):
	return settings['None']['output_path']+'counts.matrix.txt'

def gene_length_matrix_filename (settings):
	return settings['None']['output_path']+'gene.lengths.matrix.txt'

def load_settings (filename):
	"""Load the settings file
	"""
	settings = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == len(header):
			sample = row[0]
			sample_data = {}
			for el_idx, el in enumerate(header[1:]):
				sample_data[el] = row[el_idx+1]
			settings[sample] = sample_data
	return settings

def map_reads (settings, sample):
	"""Map reads using BWA
	"""
	# Make the indexes
	cmd_index = 'bwa index' + \
	            ' -p ' + bwa_index_filename(settings, sample) + \
	            ' ' + settings[sample]['fasta_file']
	print("Making index: "+cmd_index)
	subprocess.call(cmd_index, shell=True)
	# Perform the mapping
	sam_file = sam_filename(settings, sample)
	cmd_mapping = ''
	if settings[sample]['R2_fastq_file'] == '':
		cmd_mapping = 'bwa mem' + \
	                  ' ' + bwa_index_filename(settings, sample) + \
	                  ' ' + settings[sample]['R1_fastq_file'] + \
	                  ' > ' + sam_file
	else:
		cmd_mapping = 'bwa mem' + \
	                  ' ' + bwa_index_filename(settings, sample) + \
	                  ' ' + settings[sample]['R1_fastq_file'] + \
	                  ' ' + settings[sample]['R2_fastq_file'] + \
	                  ' > ' + sam_file
	print("Mapping Reads: "+cmd_mapping)
	subprocess.call(cmd_mapping, shell=True)
	# Convert to BAM for some tools
	cmd_to_bam = 'samtools view -bS' + \
	             ' ' + sam_file + \
	             ' | samtools sort' + \
	             ' -o ' + bam_filename(settings, sample, extension=True) + \
	             ' -T ' + bam_filename(settings, sample, extension=False) + \
	             ' -' + \
	             ' && samtools index ' + bam_filename(settings, sample, extension=True)
	print("Converting SAM to position sorted BAM: "+cmd_to_bam)
	subprocess.call(cmd_to_bam, shell=True)

def count_reads (settings, sample, feature='gene', attribute='Name', strand_opt='reverse'):
	"""Count reads falling in a specific feature type, group on an attribute
	"""
	# Use HTSeq to count the reads in specific features
	stranded = 'yes'
	if settings[sample]['R2_fastq_file'] == '':
		strand_opt = 'no'
	cmd_count = 'htseq-count' + \
	            ' -f bam' + \
	            ' -s ' + strand_opt + \
	            ' -a 10' + \
	            ' -m union' + \
	            ' -r pos' + \
	            ' -t ' + feature + \
	            ' -i ' + attribute + \
	            ' ' + bam_filename(settings, sample, extension=True) + \
	            ' ' + settings[sample]['gff_file'] + \
	            ' > ' + count_filename(settings, sample)
	print("Counting reads: "+cmd_count)
	subprocess.call(cmd_count, shell=True)

def mapped_reads (settings, sample):
	cmd_total = 'samtools view -c -F 4' + \
	            ' ' + bam_filename(settings, sample) + \
	            ' > ' + mapped_reads_filename(settings, sample)
	print("Total mapped reads: "+cmd_total)
	subprocess.call(cmd_total, shell=True)

def load_mapped_reads (settings, sample):
	file_in = open(mapped_reads_filename(settings, sample), 'rU')
	file_data = file_in.readlines()
	if len(file_data) > 0:
		return int(file_data[0])
	else:
		return 0

def load_gene_lengths (settings, sample):
	gene_lengths = {}
	data_reader = csv.reader(open(gene_length_filename(settings, sample), 'rU'), delimiter='\t')
	header = next(data_reader)
	for row in data_reader:
		if len(row) == 2:
			gene_lengths[row[0]] = int(row[1])
	return gene_lengths

def read_count_file (filename):
	""" Read the count file generated by HTSeq
		count_data is a dict (tag -> count)
	"""
	count_data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	for row in data_reader:
		# Check that data exists and is not reporting by HTSeq
		if len(row) == 2:
			if len(row[0])>2 and row[0][0:2] != '__':
				count_data[row[0]] = int(row[1])
	return count_data

def combine_counts (counts, sample_names):
	""" Combine a set of count dictionaries
		counts is a dictorinary of count_data where key is sample name
	"""
	full_tag_list = []
	num_of_samples = len(sample_names)
	# Generate the complete tag list (some samples won't have some tags)
	for sample in sample_names:
		full_tag_list = full_tag_list + counts[sample].keys()
	full_tag_list = list(set(full_tag_list))
	# Generate matrix zero matrix
	count_matrix = {}
	for tag in full_tag_list:
		count_matrix[tag] = [0]*num_of_samples
	# Update where count exists
	for sample_idx in range(num_of_samples):
		sample = sample_names[sample_idx]
		for tag in counts[sample].keys():
			count_matrix[tag][sample_idx] = counts[sample][tag]
	return count_matrix

def save_count_matrix (count_matrix, sample_names, filename):
	""" Save a count_matrix with the sample_names to file
	"""
	f_out = open(filename, 'w')
	# Write the header
	f_out.write( '\t'.join(['gene_name']+sample_names)+'\n' )
	for tag in sorted(count_matrix):
		count_strs = [str(x) for x in count_matrix[tag]]
		f_out.write( '\t'.join([tag]+count_strs)+'\n' )
	f_out.close()

def save_mapped_reads_matrix (mapped_reads, sample_names, filename):
	f_out = open(filename, 'w')
	f_out.write( 'sample\ttotal_mapped_reads\n' )
	for s in sample_names:
		f_out.write( s+'\t'+str(mapped_reads[s])+'\n' )
	f_out.close()

def save_gene_length_matrix (gene_lengths, filename):
	f_out = open(filename, 'w')
	f_out.write( 'gene\tlength\n' )
	seen = []
	for s in gene_lengths.keys():
		for gene in gene_lengths[s].keys():
			if gene not in seen:
				f_out.write( gene+'\t'+str(gene_lengths[s][gene])+'\n' )
				seen.append(gene)
	f_out.close()

def count_matrix (settings):
	counts = {}
	for sample in settings.keys():
		if sample != 'None':
			counts[sample] = read_count_file(count_filename(settings, sample))
	sample_names = counts.keys()
	count_matrix = combine_counts(counts, sample_names)
	save_count_matrix(count_matrix, sample_names, 
		              settings['None']['output_path']+'read_count.matrix')

def gene_lengths (settings, sample, feature='gene', attribute='Name'):
	""" Calculate the gene lengths from set of GTF references
	"""
	len_file = gene_length_filename(settings, sample)
	f_out = open(len_file, 'w')
	f_out.write('gene_name\tlength\n')
	seen = []
	data_reader = csv.reader(open(settings[sample]['gff_file'], 'rU'), delimiter='\t')
	for row in data_reader:
		if len(row) == 9 and row[2] == feature:
			attribs = row[8].split(';')
			for el in attribs:
				key = el.split('=')[0]
				value = el.split('=')[1]
				if key == attribute and value not in seen:
					gene_length = int(row[4])-int(row[3])+1
					f_out.write(value+'\t'+str(gene_length)+'\n')
					seen.append(value)
	f_out.close()

def make_profile (settings, sample):
	""" Calculate transcription profile for given regions in BED file
	"""
	if settings[sample]['R2_fastq_file'] == '':
		cmd_coverage = 'bedtools coverage -d -abam' + \
		               ' ' + bam_filename(settings, sample, extension=True) + \
		               ' -b ' + settings[sample]['bed_file'] + \
		               ' > ' + profile_filename(settings, sample)
		print("Making profile: "+cmd_coverage)
		subprocess.call(cmd_coverage, shell=True)
	else:
		cmd_fwd_coverage = 'samtools view -b -F 0x20 '+bam_filename(settings, sample, extension=True)+' |' + \
	                       ' bedtools bamtobed -i - |' + \
	                       ' bedtools coverage -d -a - -b '+settings[sample]['bed_file'] + \
	                       ' > '+profile_fwd_filename(settings, sample)
		print("Making forward profile: "+cmd_fwd_coverage)
		subprocess.call(cmd_fwd_coverage, shell=True)
		cmd_rev_coverage = 'samtools view -b -f 0x20 '+bam_filename(settings, sample, extension=True)+' |' + \
	                       ' bedtools bamtobed -i - |' + \
	                       ' bedtools coverage -d -a - -b '+settings[sample]['bed_file'] + \
	                       ' > '+profile_rev_filename(settings, sample)
		print("Making reverse profile: "+cmd_rev_coverage)
		subprocess.call(cmd_rev_coverage, shell=True)

def norm_fpkm (settings, bin_path=''):
	""" Calculate normalised RPKM/FPKM expression levels
	"""
	cmd_fpkm = 'Rscript '+bin_path+'norm_fpkm.r ' + \
	           ' ' + count_matrix_filename(settings) + \
	           ' ' + mapped_reads_matrix_filename(settings) + \
	           ' ' + gene_length_matrix_filename (settings) + \
	           ' ' + settings['None']['output_path']
	print("Generating normalised RPKM/FPKMs: "+cmd_fpkm)
	subprocess.call(cmd_fpkm, shell=True)

def de_analysis (settings, group1, group2, output_prefix, bin_path=''):
	""" Calculate DEG analysis between two groups (column numbers 1-indexed)
	"""
	cmd_deg = 'Rscript '+bin_path+'de_analysis.r ' + \
	           ' ' + count_matrix_filename (settings) + \
	           ' ' + group1 + \
	           ' ' + group2 + \
	           ' ' + mapped_reads_matrix_filename(settings) + \
	           ' ' + settings['None']['output_path'] + '/' + output_prefix
	print("Generating normalised RPKM/FPKMs: "+cmd_deg)
	subprocess.call(cmd_deg, shell=True)
