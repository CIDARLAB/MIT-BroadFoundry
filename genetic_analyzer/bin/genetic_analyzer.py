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
#   - R + edgeR package (Rscript)

# Required modules
import csv
import subprocess
import numpy as np

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
		http://seqanswers.com/forums/showthread.php?t=29399
	"""
	if settings[sample]['R2_fastq_file'] == '':
		cmd_coverage = 'bedtools coverage -d -abam' + \
		               ' ' + bam_filename(settings, sample, extension=True) + \
		               ' -b ' + settings[sample]['bed_file'] + \
		               ' > ' + profile_filename(settings, sample)
		print("Making profile: "+cmd_coverage)
		subprocess.call(cmd_coverage, shell=True)
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

def load_profiles (settings, sample):
	""" Profiles have the form of a list chr: [start_bp, end_bp, [profile_fwd],[profile_rev]]
	"""
	profiles = {}
	fwd_profile_filename = profile_fwd_filename(settings, sample)
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
				new_profile[2][int(row[3])-1] = int(row[4])
				profiles[cur_chrom].append(new_profile)
			else:
				cur_profile[0][int(row[3])-1] = int(row[4])
	rev_profile_filename = profile_rev_filename(settings, sample)
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
				cur_profile[1][int(row[3])-1] = int(row[4])
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

def characterize_promoters (settings, sample, upstream_bp=50, downstream_bp=200):
	profiles = load_profiles(settings, sample)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
						             (part_data[2]-1)-upstream_bp, part_data[3]+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
						             (part_data[2]-1)-downstream_bp, part_data[3]+upstream_bp))
				# Calculate performance
				avg_us = np.median(raw_region[0][0:upstream_bp])
				avg_ds = np.median(raw_region[0][-downstream_bp:])
				perf = avg_ds-avg_us
				char_data.append([chrom, part_name, avg_us, avg_ds, perf])
	return char_data

def characterize_terminators (settings, sample, upstream_bp=200, downstream_bp=50):
	profiles = load_profiles(settings, sample)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'terminator':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
						             (part_data[2]-1)-upstream_bp, part_data[3]+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
						             (part_data[2]-1)-downstream_bp, part_data[3]+upstream_bp))
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

def characterize_ribozymes (settings, sample, upstream_bp=25, downstream_bp=25):
	profiles = load_profiles(settings, sample)
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
						             cur_site_bp-upstream_bp, cur_site_bp+downstream_bp)
				else:
					cur_site_bp = (part_data[3])-cut_site
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
						             cur_site_bp-downstream_bp, cur_site_bp+upstream_bp))
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
							c_e = 1.0-(1.0/float(avg_ds))
						max_cut = 'Y'
					else:
						if avg_ds < avg_us:
							c_e = 0.0
						else:
							c_e = 1.0-(float(avg_us)/float(avg_ds))
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
