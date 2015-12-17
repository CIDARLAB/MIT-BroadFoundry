#!/usr/bin/env python

#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.


# Requires following software is available and in user path:
#   - BWA
#   - SAMTools
#   - HTSeq
#   - BEDTools
#   - R + edgeR package (RScript)


# Required modules
import csv
import subprocess
import numpy as np
import re
import math
import scipy.optimize


def bwa_index_filename (settings, sample):
	return settings[sample]['temp_path']+sample

def sam_filename (settings,  sample):
	return settings[sample]['temp_path']+sample+'.sam'

def bam_filename (settings,  sample, extension=True):
	bam_filename = settings[sample]['temp_path']+sample
	if extension == True:
		bam_filename += '.bam'
	return bam_filename

def fragment_dist_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fragment.distribution.txt'

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

def profile_norm_fwd_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fwd.norm.profiles.txt'

def profile_norm_rev_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.rev.norm.profiles.txt'

def count_matrix_filename (settings):
	return settings['None']['output_path']+'counts.matrix.txt'

def gene_length_matrix_filename (settings):
	return settings['None']['output_path']+'gene.lengths.matrix.txt'

def promoter_profile_perf_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.promoter.profile.perf.txt'

def terminator_profile_perf_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.terminator.profile.perf.txt'

def ribozyme_profile_perf_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.ribozyme.profile.perf.txt'

def combined_promoter_profile_perf_filename (settings):
	return settings['None']['output_path']+'promoter.profile.perf.txt'

def combined_terminator_profile_perf_filename (settings):
	return settings['None']['output_path']+'terminator.profile.perf.txt'

def combined_ribozyme_profile_perf_filename (settings):
	return settings['None']['output_path']+'ribozyme.profile.perf.txt'

def combined_fitted_promoter_perf_filename (settings, output_name):
	return settings['None']['output_path']+'fitted.promoter.perf.'+output_name+'.txt'

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

def trim_adaptors (settings, sample):
	adaptor_seq = settings[sample]['seq_adaptor']
	cmd_index = 'bowtie-build index' + \
				' -p ' + bwa_index_filename(settings, sample) + \
				' ' + settings[sample]['fasta_file']
	print("Making index: "+cmd_index)
	subprocess.call(cmd_index, shell=True)


def map_ribo_reads (settings, sample):
	cmd_index = 'bowtie-build index' + \
				' -p ' + bwa_index_filename(settings, sample) + \
				' ' + settings[sample]['fasta_file']
	print("Making index: "+cmd_index)
	subprocess.call(cmd_index, shell=True)

	cmd_index = 'bowtie ' + \
				' -p ' + bwa_index_filename(settings, sample) + \
				' ' + settings[sample]['fasta_file']
	print("Perform mapping: "+cmd_index)
	subprocess.call(cmd_index, shell=True)


def map_reads (settings, sample):
	"""Map reads using BWA-MEM
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
	print("Mapping Reads (BWM-MEM): "+cmd_mapping)
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
	if settings[sample]['R2_fastq_file'] == '' and strand_opt == 'reverse':
		strand_opt = 'yes'
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
		http://seqanswers.com/forums/showthread.php?t=29399
	"""
	if settings[sample]['R2_fastq_file'] == '':
		# https://www.biostars.org/p/14378/
		fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
		cmd_fwd_coverage = 'samtools view -b -F 20 '+bam_filename(settings, sample, extension=True)+' > '+fwd_filename+' && '+\
			'samtools index '+fwd_filename+' && '+\
			'bedtools coverage -d -abam '+fwd_filename+' -b '+settings[sample]['bed_file'] + \
			' > '+profile_fwd_filename(settings, sample)
		print("Making forward profile: "+cmd_fwd_coverage)
		subprocess.call(cmd_fwd_coverage, shell=True)
		rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
		cmd_rev_coverage = 'samtools view -b -f 16 '+bam_filename(settings, sample, extension=True)+' > '+rev_filename+' && '+\
			'samtools index '+rev_filename+' && '+\
			'bedtools coverage -d -abam '+rev_filename+' -b '+settings[sample]['bed_file'] + \
			' > '+profile_rev_filename(settings, sample)
		print("Making reverse profile: "+cmd_rev_coverage)
		subprocess.call(cmd_rev_coverage, shell=True)

	else:
		fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
		fwd1_filename = bam_filename(settings, sample, extension=False) + '.fwd.1.bam'
		fwd2_filename = bam_filename(settings, sample, extension=False) + '.fwd.2.bam'
		cmd_fwd_coverage = 'samtools view -b -f 83 '+bam_filename(settings, sample, extension=True)+' > '+fwd1_filename+' && '+\
			'samtools index '+fwd1_filename+' && '+\
			'samtools view -b -f 163 '+bam_filename(settings, sample, extension=True)+' > '+fwd2_filename+' && '+\
			'samtools index '+fwd2_filename+' && '+\
			'samtools merge -f '+fwd_filename+' '+fwd1_filename+' '+fwd2_filename+' && '+\
			'samtools index '+fwd_filename+' && '+\
			'bedtools coverage -d -abam '+fwd_filename+' -b '+settings[sample]['bed_file'] + \
			' > '+profile_fwd_filename(settings, sample)
		print("Making forward profile: "+cmd_fwd_coverage)
		subprocess.call(cmd_fwd_coverage, shell=True)
		rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
		rev1_filename = bam_filename(settings, sample, extension=False) + '.rev.1.bam'
		rev2_filename = bam_filename(settings, sample, extension=False) + '.rev.2.bam'
		cmd_rev_coverage = 'samtools view -b -f 99 '+bam_filename(settings, sample, extension=True)+' > '+rev1_filename+' && '+\
			'samtools index '+rev1_filename+' && '+\
			'samtools view -b -f 147 '+bam_filename(settings, sample, extension=True)+' > '+rev2_filename+' && '+\
			'samtools index '+rev2_filename+' && '+\
			'samtools merge -f '+rev_filename+' '+rev1_filename+' '+rev2_filename+' && '+\
			'samtools index '+rev_filename+' && '+\
			'bedtools coverage -d -abam '+rev_filename+' -b '+settings[sample]['bed_file'] + \
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
			   ' ' + settings['None']['output_path'] + output_prefix
	print("Generating normalised RPKM/FPKMs: "+cmd_deg)
	subprocess.call(cmd_deg, shell=True)

########## CHARACTERIZATION ##########

def load_gff (settings, sample):
	gff = {}
	data_reader = csv.reader(open(settings[sample]['gff_file'], 'rU'), delimiter='\t')
	# Process each line
	for row in data_reader:
		if len(row) == 9:
			chromo = row[0]
			part_type = row[2]
			start_bp = int(row[3])
			end_bp = int(row[4])
			part_dir = row[6]
			part_attribs = {}
			split_attribs = row[8].split(';')
			part_name = None
			for attrib in split_attribs:
				key_value = attrib.split('=')
				if len(key_value) == 2:
					if key_value[0] == 'Name':
						part_name = key_value[1]
					else:
						part_attribs[key_value[0]] = key_value[1]
			if part_name != None:
				if chromo not in gff.keys():
					gff[chromo] = {}
				gff[chromo][part_name] = [part_type, part_dir, start_bp, end_bp, part_attribs]
	return gff

def load_profiles (settings, sample, normed=False):
	""" Profiles have the form of a list chr: [start_bp, end_bp, [profile_fwd],[profile_rev]]
	"""
	profiles = {}
	fwd_profile_filename = profile_fwd_filename(settings, sample)
	if normed == True:
		fwd_profile_filename = profile_norm_fwd_filename(settings, sample)
	data_reader = csv.reader(open(fwd_profile_filename, 'rU'), delimiter='\t')
	# Process each line in fwd profile
	for row in data_reader:
		if len(row) == 5:
			cur_chrom = row[0]
			if cur_chrom not in profiles.keys():
				profiles[cur_chrom] = []
			cur_start_bp = int(row[1])
			cur_end_bp = int(row[2])
			cur_profile = find_profile(profiles, cur_chrom, cur_start_bp, cur_end_bp)
			if cur_profile == None:
				new_profile = [cur_start_bp, cur_end_bp, np.zeros(cur_end_bp-cur_start_bp), np.zeros(cur_end_bp-cur_start_bp)]
				new_profile[2][int(row[3])-1] = float(row[4])
				profiles[cur_chrom].append(new_profile)
			else:
				cur_profile[0][int(row[3])-1] = float(row[4])
	rev_profile_filename = profile_rev_filename(settings, sample)
	if normed == True:
		rev_profile_filename = profile_norm_rev_filename(settings, sample)
	data_reader = csv.reader(open(rev_profile_filename, 'rU'), delimiter='\t')
	# Process each line in fwd profile
	for row in data_reader:
		if len(row) == 5:
			cur_chrom = row[0]
			if cur_chrom not in profiles.keys():
				profiles[cur_chrom] = []
			cur_start_bp = int(row[1])
			cur_end_bp = int(row[2])
			cur_profile = find_profile(profiles, cur_chrom, cur_start_bp, cur_end_bp)
			if cur_profile != None:
				cur_profile[1][int(row[3])-1] = float(row[4])
	return profiles

def find_profile (profiles, chrom, start_bp, end_bp):
	if chrom in profiles.keys():
		for el in profiles[chrom]:
			if el[0] == start_bp and el[1] == end_bp:
				return [el[2], el[3]]
	return None

def extract_profile_region (profiles, chrom, start_bp, end_bp):
	region = None
	if chrom in profiles.keys():
		for profile in profiles[chrom]:
			full_chrom = False
			if profile[0] == 0 and profile[1] == len(profile[2]):
				full_chrom = True
			if full_chrom == True:
				fwd_profile = list(profile[2])
				rev_profile = list(profile[3])
				profile_len = len(fwd_profile)
				ext_start_fwd = []
				ext_end_fwd = []
				ext_start_rev = []
				ext_end_rev = []
				# The region will exist
				if start_bp < 0:
					# extend the profile at start
					ext_start_fwd = fwd_profile[start_bp:]
					ext_start_rev = rev_profile[start_bp:]
				if end_bp > profile_len:
					# extend the profile at end
					ext_end_fwd = fwd_profile[:(end_bp-profile_len)]
					ext_end_rev = rev_profile[:(end_bp-profile_len)]
				new_start_bp = start_bp
				new_end_bp = end_bp
				if ext_start_fwd != []:
					new_start_bp = 0
					new_end_bp = end_bp+len(ext_start_fwd)
				new_fwd_profile = ext_start_fwd+fwd_profile+ext_end_fwd
				new_rev_profile = ext_start_rev+rev_profile+ext_end_rev
				region = [new_fwd_profile[new_start_bp:new_end_bp], 
						  new_rev_profile[new_start_bp:new_end_bp]]
				break
	return region

def reverse_region (region):
	return [region[1][::-1], region[0][::-1]]

def characterize_promoter_units (settings, sample, upstream_bp=25, downstream_skip_bp=0, downstream_bp=50, normed=False):
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter_unit':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-upstream_bp, part_data[3]+downstream_skip_bp+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-downstream_skip_bp-downstream_bp, part_data[3]+upstream_bp))
				# Calculate performance
				avg_us = np.median(raw_region[0][0:upstream_bp])
				avg_ds = np.median(raw_region[0][-downstream_bp:])
				perf = avg_ds-avg_us
				char_data.append([chrom, part_name, avg_us, avg_ds, perf])
	return char_data

def characterize_promoters (settings, sample, upstream_bp=25, downstream_skip_bp=0, downstream_bp=50, normed=False):
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-upstream_bp, part_data[3]+downstream_skip_bp+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-downstream_skip_bp-downstream_bp, part_data[3]+upstream_bp))
				# Calculate performance
				avg_us = np.median(raw_region[0][0:upstream_bp])
				avg_ds = np.median(raw_region[0][-downstream_bp:])
				perf = avg_ds-avg_us
				char_data.append([chrom, part_name, avg_us, avg_ds, perf])
	return char_data

def characterize_terminators (settings, sample, upstream_bp=50, upstream_skip_bp=0, downstream_bp=25, normed=False):
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'terminator':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-upstream_skip_bp-upstream_bp, part_data[3]+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-downstream_bp, part_data[3]+upstream_skip_bp+upstream_bp))
				# Calculate performance
				avg_us = np.mean(raw_region[0][0:upstream_bp])
				avg_ds = np.mean(raw_region[0][-downstream_bp:])
				max_term = 'N'
				t_e = 0.0
				if avg_us == 0.0:
					max_term = 'Y'
					t_e = 0.0
				else:
					if avg_ds == 0:
						t_e = 1.0/float(avg_us)
						max_term = 'Y'
					else:
						t_e = float(avg_ds)/float(avg_us)
				if t_e != 0.0:
					t_s = 1.0/t_e
				else:
					t_s = -1.0
				char_data.append([chrom, part_name, avg_us, avg_ds, t_e, t_s, max_term])
	return char_data

def characterize_ribozymes (settings, sample, upstream_bp=50, downstream_skip_bp=0, downstream_bp=50, normed=False):
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'ribozyme':
				cut_site = 0
				if 'cut_site' in part_data[4].keys():
					cut_site = int(part_data[4]['cut_site'])
				if part_data[1] == '+':
					cur_site_bp = (part_data[2]-1)+cut_site
					raw_region = extract_profile_region(profiles, chrom, 
									 cur_site_bp-upstream_bp, cur_site_bp+downstream_skip_bp+downstream_bp)
				else:
					cur_site_bp = (part_data[3])-cut_site
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 cur_site_bp-downstream_skip_bp-downstream_bp, cur_site_bp+upstream_bp))
				# Calculate performance
				avg_us = np.mean(raw_region[0][0:upstream_bp])
				avg_ds = np.mean(raw_region[0][-downstream_bp:])
				max_cut = 'N'
				c_e = 0.0
				if avg_ds == 0.0:
					c_e = 0.0
					max_cut = 'Y'
				else:
					if avg_us == 0:
						if avg_ds < avg_us:
							c_e = 0.0
						else:
							#c_e = 1.0-(1.0/float(avg_ds))
							c_e = ((float(avg_ds)+1.0)/1.0) - 1.0
						max_cut = 'Y'
					else:
						if avg_ds < avg_us:
							c_e = 0.0
						else:
							#c_e = 1.0-(float(avg_us)/float(avg_ds))
							c_e = (float(avg_ds)/float(avg_us)) - 1.0
				char_data.append([chrom, part_name, avg_us, avg_ds, c_e, max_cut, cut_site])
	return char_data

def save_characterization_data (settings, sample, data, part_type=None):
	if part_type == 'promoter':
		filename = promoter_profile_perf_filename(settings, sample)
		f_out = open(filename, 'w')
		f_out.write( 'sample\tchromosome\tpart_name\treads_us\treads_ds\treads_strength\n' )
		for d in data:
			f_out.write( sample+'\t'+'\t'.join([str(x) for x in d])+'\n' )
		f_out.close()
	if part_type == 'terminator':
		filename = terminator_profile_perf_filename(settings, sample)
		f_out = open(filename, 'w')
		f_out.write( 'sample\tchromosome\tpart_name\treads_us\treads_ds\tt_e\tt_s\tmax_term\n' )
		for d in data:
			f_out.write( sample+'\t'+'\t'.join([str(x) for x in d])+'\n' )
		f_out.close()
	if part_type == 'ribozyme':
		filename = ribozyme_profile_perf_filename(settings, sample)
		f_out = open(filename, 'w')
		f_out.write( 'sample\tchromosome\tpart_name\treads_us\treads_ds\tc_e\tmax_cut\tcut_site\n' )
		for d in data:
			f_out.write( sample+'\t'+'\t'.join([str(x) for x in d])+'\n' )
		f_out.close()

def combine_promoter_characterizations (settings, samples):
	data = {}
	for s in samples:
		filename = promoter_profile_perf_filename(settings, s)
		data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
		header = next(data_reader)
		for row in data_reader:
			if row[1] not in data.keys():
				data[row[1]] = {}
			chrom_data = data[row[1]]
			if row[2] not in chrom_data.keys():
				chrom_data[row[2]] = []
			chrom_part_data = chrom_data[row[2]]
			chrom_part_data.append([row[0]]+row[3:])
	f_out = open(combined_promoter_profile_perf_filename(settings), 'w')
	f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\treads_strength\n')
	for chrom in sorted(data.keys()):
		chrom_data = data[chrom]
		for part in sorted(chrom_data.keys()):
			chrom_part_data = chrom_data[part]
			for data_rec in chrom_part_data:
				f_out.write( '\t'.join([chrom, part]+data_rec)+'\n' )
	f_out.close()

def combine_terminator_characterizations (settings, samples):
	data = {}
	for s in samples:
		filename = terminator_profile_perf_filename(settings, s)
		data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
		header = next(data_reader)
		for row in data_reader:
			if row[1] not in data.keys():
				data[row[1]] = {}
			chrom_data = data[row[1]]
			if row[2] not in chrom_data.keys():
				chrom_data[row[2]] = []
			chrom_part_data = chrom_data[row[2]]
			chrom_part_data.append([row[0]]+row[3:])
	f_out = open(combined_terminator_profile_perf_filename(settings), 'w')
	f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\tt_e\tt_s\tmax_term\n')
	for chrom in sorted(data.keys()):
		chrom_data = data[chrom]
		for part in sorted(chrom_data.keys()):
			chrom_part_data = chrom_data[part]
			for data_rec in chrom_part_data:
				f_out.write( '\t'.join([chrom, part]+data_rec)+'\n' )
	f_out.close()

def combine_ribozyme_characterizations (settings, samples):
	data = {}
	for s in samples:
		filename = ribozyme_profile_perf_filename(settings, s)
		data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
		header = next(data_reader)
		for row in data_reader:
			if row[1] not in data.keys():
				data[row[1]] = {}
			chrom_data = data[row[1]]
			if row[2] not in chrom_data.keys():
				chrom_data[row[2]] = []
			chrom_part_data = chrom_data[row[2]]
			chrom_part_data.append([row[0]]+row[3:])
	f_out = open(combined_ribozyme_profile_perf_filename(settings), 'w')
	f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\tc_e\tmax_cut\tcut_site\n')
	for chrom in sorted(data.keys()):
		chrom_data = data[chrom]
		for part in sorted(chrom_data.keys()):
			chrom_part_data = chrom_data[part]
			for data_rec in chrom_part_data:
				f_out.write( '\t'.join([chrom, part]+data_rec)+'\n' )
	f_out.close()

########## FRAGMENT DISTRIBUTIONS ##########

def fragment_length_dists (settings, sample, reads_to_sample=1000000):
	""" Adapted from get_insert_size.py (Wei Li)
	"""
	frag_file = fragment_dist_filename(settings, sample)
	sam_file = sam_filename(settings, sample)
	f_in = open(sam_file, 'rU')
	plrdlen={}
	plrdspan={}
	objmrl=re.compile('([0-9]+)M$')
	objmtj=re.compile('NH:i:(\d+)')
	nline=0
	with open(sam_file, 'rU') as ins:
		for lines in ins:
			field=lines.strip().split()
			nline=nline+1
			if nline >= reads_to_sample:
				break
			if len(field)<12:
				continue
			try:
				mrl=objmrl.match(field[5])
				if mrl==None: # ignore non-perfect reads
					continue
				readlen=int(mrl.group(1))
				if readlen in plrdlen.keys():
					plrdlen[readlen]=plrdlen[readlen]+1
				else:
					plrdlen[readlen]=1
				if field[6]!='=':
					continue
				dist=int(field[8])
				if dist<=0: # ignore neg dist
					continue
				mtj=objmtj.search(lines)
				if dist in plrdspan.keys():
					plrdspan[dist]=plrdspan[dist]+1
				else:
					plrdspan[dist]=1
			except ValueError:
				continue
	f_out = open(frag_file, 'w')
	for k in sorted(plrdspan.keys()):
		f_out.write(str(k)+'\t'+str(plrdspan[k])+'\n')

def load_norm_factor (settings, sample):
	norm_facs = {}
	norm_fac_file = settings['None']['output_path']+'norm.factors.matrix.txt'
	data_reader = csv.reader(open(norm_fac_file, 'rU'), delimiter='\t')
	# Ignore the header
	data_reader.next()
	# Process each line
	for row in data_reader:
		if len(row) == 3:
			norm_facs[row[0]] = [float(row[1]), float(row[2])]
	return (norm_facs[sample][0]*norm_facs[sample][1])/1000000.0

def load_fragmentation_dist (settings, sample, max_frag_len=1000):
	frag_dist = np.zeros(max_frag_len+1)
	frag_file = fragment_dist_filename(settings, sample)
	data_reader = csv.reader(open(frag_file, 'rU'), delimiter='\t')
	# Process each line
	for row in data_reader:
		frag_len = int(row[0])
		frag_count = int(row[1])
		if frag_len <= max_frag_len:
			frag_dist[frag_len] = frag_count
	return frag_dist

def correction_prob (frag_dist, transcript_len, bp_idx):
	frag_sum = 0.0
	for l in range(len(frag_dist)):
		frag_sum += frag_dist[l] * l # min([l, transcript_len-l+1])
	frag_N = 0.0
	for l in range(len(frag_dist)):
		frag_N += frag_dist[l] * min([l, bp_idx, 
			                          (transcript_len-bp_idx+1),
			                          (transcript_len-l+1)])
	return frag_N/frag_sum

def correction_profile (frag_dist, transcript_len):
	corr_profile = np.ones(transcript_len)
	for idx in range(transcript_len):
		corr_prob = correction_prob(frag_dist, transcript_len, idx+1)
		if corr_prob > 0.0:
			corr_profile[idx] = 1.0/corr_prob
	return corr_profile

def normalise_profiles (settings, sample, correction=False, baseline_us_bp=5):
	""" Will normalise the profiles using TMM factors and correct for edge effects using fragment distribution
	"""
	# Load the profiles for the sample
	profiles = load_profiles(settings, sample)
	# If correction of edge effect required
	if correction == True:
		frag_dist = load_fragmentation_dist(settings, sample)
		gff = load_gff (settings, sample)
		for chrom in gff.keys():
			for part_name in gff[chrom].keys():
				part_data = gff[chrom][part_name]
				if part_data[0] == 'transcript':
					trans_start = part_data[2]-1
					trans_end = part_data[3]
					# Update each profile that contains the transcript
					for p_data_idx in range(len(profiles[chrom])):
						p_data = profiles[chrom][p_data_idx]
						if p_data[0] <= trans_start and p_data[1] >= trans_end:
							trans_len = trans_end-trans_start

							# Calc corr_prob_dist (symetrical)
							corr_prob_dist = correction_profile(frag_dist, trans_len)

							if part_data[1] == '+':
								rel_start = trans_start-p_data[0]
								rel_end = trans_end-p_data[0]
								# Calc baseline to subtract
								bl_start = rel_start-baseline_us_bp
								baseline_reads = 0.0
								if bl_start >= 0:
									baseline_reads = np.mean(p_data[2][bl_start:rel_start])
								# Calculate profile corrected for edge effects and update
								corr_profile = corr_prob_dist*(p_data[3][rel_start:rel_end]-baseline_reads).clip(min=0.0)
								profiles[chrom][p_data_idx][2][rel_start:rel_end] = corr_profile+baseline_reads
							else:
								rel_start = trans_start-p_data[0]
								rel_end = trans_end-p_data[0]
								# Calc baseline to subtract
								bl_start = rel_end+baseline_us_bp
								baseline_reads = 0.0
								if bl_start >= p_data[1]-p_data[0]:
									baseline_reads = np.mean(p_data[3][rel_end:bl_start])
								# Calculate profile corrected for edge effects and update
								corr_profile = corr_prob_dist*(p_data[3][rel_start:rel_end]-baseline_reads).clip(min=0.0)
								profiles[chrom][p_data_idx][3][rel_start:rel_end] = corr_profile+baseline_reads

	# Save normalised profiles to file
	norm_fac = load_norm_factor(settings, sample)
	fwd_file = profile_norm_fwd_filename(settings, sample)
	rev_file = profile_norm_rev_filename(settings, sample)
	f_out_fwd = open(fwd_file, 'w')
	f_out_rev = open(rev_file, 'w')
	for chrom in profiles.keys():
		for region in profiles[chrom]:
			cur_start_bp = region[0]
			cur_end_bp = region[1]
			cur_fwd_profile = region[2]
			cur_rev_profile = region[3]
			for idx in range(len(cur_fwd_profile)):
				cur_fwd_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx+1), str(cur_fwd_profile[idx]/norm_fac)]
				cur_rev_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx+1), str(cur_rev_profile[idx]/norm_fac)]
				f_out_fwd.write('\t'.join(cur_fwd_data)+'\n')
				f_out_rev.write('\t'.join(cur_rev_data)+'\n')
	f_out_fwd.close()
	f_out_rev.close()

########## FITTING OF RESPONSE FUNCTIONS ##########

def load_promoter_characterizations (filename, samples):
	pro_data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	data_reader.next()
	# Process each line
	for row in data_reader:
		cur_chrom = row[0]
		cur_pu = row[1]
		cur_sample = row[2]
		cur_perf = float(row[5])
		if cur_sample in samples:
			if cur_sample not in pro_data.keys():
				pro_data[cur_sample] = {}
			if cur_chrom not in pro_data[cur_sample].keys():
				pro_data[cur_sample][cur_chrom] = {}
			pro_data[cur_sample][cur_chrom][cur_pu] = cur_perf
	return pro_data


def hill_func (x, Pmin, Pmin_inc, K, n, repress=False):
	if x < 0.0:
		x = 0.0
	if repress == True:
		return Pmin + (Pmin+Pmin_inc)*( math.pow(K,n) / (math.pow(K,n)+math.pow(x,n)) )
	else: 
		return Pmin + (Pmin+Pmin_inc)*( math.pow(x,n) / (math.pow(K,n)+math.pow(x,n)) )


def extract_fit_params (x, P_names, P_types):
	# Extract the parameters into more user friendly dict
	fit_params = {}
	cur_idx = 0
	for p_idx in range(len(P_names)):
		p_name = P_names[p_idx]
		p_type = P_types[p_idx]
		if p_type == 'induced':
			# If an induced promoter then 2 parameters
			fit_params[p_name] = {}
			fit_params[p_name]['Pmin'] = x[cur_idx]
			fit_params[p_name]['Pmin_inc'] = x[cur_idx+1]
			cur_idx += 2
		else:
			# If a repessed promoter then 4 parameters
			fit_params[p_name] = {}
			fit_params[p_name]['Pmin'] = x[cur_idx]
			fit_params[p_name]['Pmin_inc'] = x[cur_idx+1]
			fit_params[p_name]['K'] = x[cur_idx+2]
			fit_params[p_name]['n'] = x[cur_idx+3]
			cur_idx += 4
	return fit_params


def promoter_unit_err_func (x, exp_data, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names, P_types, fac_reduce_size):		
	# Extract the parameters into more user friendly dict
	fit_params = extract_fit_params(x, P_names, P_types)
	# Only calculate error based on samples specified
	err_diffs = []
	for sample in exp_data.keys():
		sample_data = exp_data[sample]
		exp_val = exp_data[sample][chrom_to_fit][PU_to_fit]
		# Calculate fitted output assuming additive fluxes
		fit_outs = []
		cur_param_idx = 0
		for p_idx in range(len(P_names)):
			cur_p = P_names[p_idx]
			# Check the type of promoter and calculate 
			if P_types[p_idx] == 'induced':
				# Inducible promoter
				input_val = extract_key_vals(PU_inputs[p_idx])[sample]
				if input_val == 1.0:
					fit_outs.append(fit_params[cur_p]['Pmin']+fit_params[cur_p]['Pmin_inc'])
				else:
					fit_outs.append(fit_params[cur_p]['Pmin'])
			else:
				# Repressor promoter
				input_val = exp_data[sample][chrom_inputs[p_idx]][PU_inputs[p_idx]]
				fit_outs.append( hill_func(input_val, 
					                       fit_params[cur_p]['Pmin'], 
					                       fit_params[cur_p]['Pmin_inc'],
					                       fit_params[cur_p]['K'],
					                       fit_params[cur_p]['n'],
					                       repress=True) )
		fit_val = np.sum(fit_outs)
		# Numbers are generally in read-depths (reduce to reduce squared values being too large)
		err_diffs.append((exp_val-fit_val)/fac_reduce_size)
	# Return SSE
	return np.sum(np.power(err_diffs, 2.0))


def extract_key_vals (data_str):
	split_str = data_str.split(':')
 	key_vals = {}
	for el in split_str:
		key_val_split = el.split('>')
		if len(key_val_split) == 2:
			key_vals[key_val_split[0]] = float(key_val_split[1])
	return key_vals


def extract_list (data_str):
	split_str = data_str.split(',')
	return split_str
	

def fit_promoter_response_functions (settings, samples, output_name, fac_reduce_size=1000.0):
	# All the data needed for the fitting
	promoter_filename = combined_promoter_profile_perf_filename(settings)
	pro_data = load_promoter_characterizations(promoter_filename, samples)
	gff = load_gff(settings, samples[0])
	# Somewhere to save the results
	fitted_pro_params = {}
	# Parameters required for each fiting
	chrom_to_fit = ''
	PU_to_fit = ''
	chrom_inputs = []
	PU_inputs = []
	P_names = []
	P_types = []
	# Cycle through each promoter unit and fit the internal promoter functions
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter_unit':
				part_attribs = part_data[-1]
				# Populate all the parameters of the promoter unit
				chrom_to_fit = chrom
				PU_to_fit = part_name
				chrom_inputs = extract_list(part_attribs['chrom_inputs'])
				P_names = extract_list(part_attribs['promoter_names'])
				P_types = extract_list(part_attribs['promoter_types'])
				PU_inputs = extract_list(part_attribs['promoter_unit_inputs'])
				# Calculate the number of parameters we have
				num_of_params = 0
				for p_idx in range(len(P_names)):
					if P_types[p_idx] == 'induced':
						num_of_params += 2
					else:
						num_of_params += 4
				# Parameters for fit
				x0 = np.zeros(num_of_params)
				# Some contraints (keep everything positive)
				bnds = []
				cur_param_idx = 0
				for p_idx in range(len(P_names)):
					if P_types[p_idx] == 'induced':
						# Initial conditions (start realistic to improve fitting)
						x0[cur_param_idx] = 0.0
						x0[cur_param_idx+1] = 1000.0
						# Set bounds
						bnds.append((0.0, None))
						bnds.append((0.0, None))
						cur_param_idx += 2
					else:
						# Initial conditions (start realistic to improve fitting)
						x0[cur_param_idx] = 0.0
						x0[cur_param_idx+1] = 1000.0
						x0[cur_param_idx+2] = 100.0
						x0[cur_param_idx+3] = 2.0
						# Set bounds
						bnds.append((0.0, None))
						bnds.append((0.0, None))
						bnds.append((1.0, None))
						bnds.append((1.0, 3.9))
						cur_param_idx += 4
				# methods = BFGS, nelder-mead, Powell, TNC, SLSQP,  L-BFGS-B
				res = scipy.optimize.minimize(promoter_unit_err_func, x0, args=(pro_data, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names, P_types, fac_reduce_size),
											  bounds=bnds,
					                          method='SLSQP', jac=False,
					                          options={'disp': True, 'maxiter': 10000})
				# Save the fitted parameters
				cur_param_idx = 0
				for p_idx in range(len(P_names)):
					p_name = P_names[p_idx]
					if P_types[p_idx] == 'induced':
						fitted_pro_params[p_name] = {}
						fitted_pro_params[p_name]['Pmin'] = res.x[cur_param_idx]
						fitted_pro_params[p_name]['Pmin_inc'] = res.x[cur_param_idx+1]
						cur_param_idx += 2
					else:
						fitted_pro_params[p_name] = {}
						fitted_pro_params[p_name]['Pmin'] = res.x[cur_param_idx]
						fitted_pro_params[p_name]['Pmin_inc'] = res.x[cur_param_idx+1]
						fitted_pro_params[p_name]['K'] = res.x[cur_param_idx+2]
						fitted_pro_params[p_name]['n'] = res.x[cur_param_idx+3]
						cur_param_idx += 4
	# Save the results to file	
	out_filename = combined_fitted_promoter_perf_filename(settings, output_name)
	f_out = open(out_filename, 'w')
	for p_name in fitted_pro_params.keys():
		f_out.write(p_name)
		for param in fitted_pro_params[p_name].keys():
			f_out.write('\t'+param+'\t'+str(fitted_pro_params[p_name][param]))
		f_out.write('\n')
	f_out.close()
