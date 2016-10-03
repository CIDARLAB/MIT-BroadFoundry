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
#   - GATK3
#   - Picard-tools
#   - Java


# Required modules
import sys
import imp
import csv
import subprocess
import re
import math
import os
from collections import OrderedDict
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.optimize
import numpy as np


def bwa_index_filename(settings, sample):
    return settings[sample]['temp_path'] + sample


def sam_filename(settings, sample):
    return settings[sample]['temp_path'] + sample + '.sam'


def bam_filename(settings, sample, extension=True):
    bam_filename = settings[sample]['temp_path'] + sample
    if extension == True:
        bam_filename += '.bam'
    return bam_filename


def bam_duplicate_name(settings, sample, extension=True):
    bam_duplicate_name = settings[sample]['temp_path'] + sample + '.mdup'
    if extension == True:
        bam_duplicate_name += '.bam'
    return bam_duplicate_name


def fragment_dist_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.fragment.distribution.txt'


def count_filename(settings, sample, feature):
    return settings[sample]['output_path'] + sample + '.' + feature + '.counts.txt'


def mapped_reads_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.mapped.reads.txt'


def mapped_reads_matrix_filename(settings):
    return settings['None']['output_path'] + 'mapped.reads.matrix.txt'


def gene_length_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.gene.lengths.txt'


def profile_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.profiles.txt'


def profile_fwd_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.fwd.profiles.txt'


def profile_rev_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.rev.profiles.txt'


def profile_norm_fwd_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.fwd.norm.profiles.txt'


def profile_norm_rev_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.rev.norm.profiles.txt'


def count_matrix_filename(settings):
    return settings['None']['output_path'] + 'counts.matrix.txt'


def normed_counts_filename(settings):
    return settings['None']['output_path'] + 'fpkm.normed.matrix.txt'


def gene_length_matrix_filename(settings):
    return settings['None']['output_path'] + 'gene.lengths.matrix.txt'


def promoter_profile_perf_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.promoter.profile.perf.txt'


def terminator_profile_perf_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.terminator.profile.perf.txt'


def ribozyme_profile_perf_filename(settings, sample):
    return settings[sample]['output_path'] + sample + '.ribozyme.profile.perf.txt'


def combined_promoter_profile_perf_filename(settings):
    return settings['None']['output_path'] + 'promoter.profile.perf.txt'


def combined_terminator_profile_perf_filename(settings):
    return settings['None']['output_path'] + 'terminator.profile.perf.txt'


def combined_ribozyme_profile_perf_filename(settings):
    return settings['None']['output_path'] + 'ribozyme.profile.perf.txt'


def combined_fitted_promoter_perf_filename(settings, output_name):
    return settings['None']['output_path'] + 'fitted.promoter.perf.' + output_name + '.txt'


def promoter_reu_filename(settings, sample):
    """Load promoter reu filename
    Modify /rna-seq/ if the destination of RNA seq related project specific REU data changes.

    Args:
        settings: (dict) Described in load_settings()
        sample: string, which sample to look for within settings dict

    Returns:
        (str) File path of location for promoter reu data.
    """
    synthetic_construct = settings[sample]['fasta_file'].split('/')[-1].split('.')[0]
    prefix = '/'.join(settings[sample]['fasta_file'].split('/')[:4])
    return str(prefix) + '/rna-seq/' +  str(synthetic_construct) + '/' + 'promoter_reu.txt'


def terminator_reu_filename(settings, sample):
    """Load terminator reu filename
    Modify /rna-seq/ if the destination of RNA seq related project specific REU data changes.

    Args:
        settings: (dict) Described in load_settings()
        sample: string, which sample to look for within settings dict

    Returns:
        (str) File path of location for terminator reu data.
    """
    synthetic_construct = settings[sample]['fasta_file'].split('/')[-1].split('.')[0]
    prefix = '/'.join(settings[sample]['fasta_file'].split('/')[:4])
    return str(prefix) + '/rna-seq/' +  str(synthetic_construct) + '/' + 'terminator_reu.txt'


def vcf_name(settings, sample):
    return settings[sample]['output_path'] + sample + '.vcf'


def context_filename(settings):
    return settings['None']['output_path'] + 'context_data.txt'


def bam_readgroup_name(settings, sample, extension=True):
    bam_readgroup_name = settings[sample]['temp_path'] + sample + '.mdup.rg'
    if extension == True:
        bam_readgroup_name += '.bam'
    return bam_readgroup_name


def load_settings(filename):
    """Load the settings file
    """
    settings = {}
    try:
	    data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
    except IOError:
        return sys.exit(1)
    # Ignore header
    header = next(data_reader)
    # Process each line
    for row in data_reader:
        if len(row) == len(header):
            sample = row[0]
            sample_data = {}
            for el_idx, el in enumerate(header[1:]):
                sample_data[el] = row[el_idx + 1]
            settings[sample] = sample_data
    return settings


def load_features(settings, sample):
    gff = load_gff(settings, sample)
    features = set()
    for chrom, part in gff.items():
        for g in gff[chrom].values():
            features.add(g[0])
    return features

def trim_adaptors(settings, sample):
    adaptor_seq = settings[sample]['seq_adaptor']
    cmd_index = 'bowtie-build index' + \
                ' -p ' + bwa_index_filename(settings, sample) + \
                ' ' + settings[sample]['fasta_file']
    print("Making index: " + cmd_index)
    subprocess.call(cmd_index, shell=True)


def map_ribo_reads(settings, sample):
    cmd_index = 'bowtie-build index' + \
                ' -p ' + bwa_index_filename(settings, sample) + \
                ' ' + settings[sample]['fasta_file']
    print("Making index: " + cmd_index)
    subprocess.call(cmd_index, shell=True)

    cmd_index = 'bowtie ' + \
                ' -p ' + bwa_index_filename(settings, sample) + \
                ' ' + settings[sample]['fasta_file']
    print("Perform mapping: " + cmd_index)
    subprocess.call(cmd_index, shell=True)


def map_reads(settings, sample):
    """Map reads using BWA-MEM
    """
    # Make the indexes
    cmd_index = 'bwa index' + \
                ' -p ' + bwa_index_filename(settings, sample) + \
                ' ' + settings[sample]['fasta_file']
    print("Making index: " + cmd_index)
    status3 = subprocess.call(cmd_index, shell=True)
    if status3 != 0:
        return 1
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
    print("Mapping Reads (BWM-MEM): " + cmd_mapping)
    status1 = subprocess.call(cmd_mapping, shell=True)
    if status1 != 0:
        return 1
    # Convert to BAM for some tools
    cmd_to_bam = 'samtools view -bS' + \
                 ' ' + sam_file + \
                 ' | samtools sort' + \
                 ' -o ' + bam_filename(settings, sample, extension=True) + \
                 ' -T ' + bam_filename(settings, sample, extension=False) + \
                 ' -' + \
                 ' && samtools index ' + bam_filename(settings, sample, extension=True)
    print("Converting SAM to position sorted BAM: " + cmd_to_bam)
    status2 = subprocess.call(cmd_to_bam, shell=True)
    return status2


def count_reads(settings, sample, feature='gene', attribute='Name', strand_opt='reverse'):
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
                ' > ' + count_filename(settings, sample, feature)
    print("Counting reads: " + cmd_count)
    status = subprocess.call(cmd_count, shell=True)
    return status


def mapped_reads(settings, sample):
    cmd_total = 'samtools view -c -F 4' + \
                ' ' + bam_filename(settings, sample) + \
                ' > ' + mapped_reads_filename(settings, sample)
    print("Total mapped reads: " + cmd_total)
    status = subprocess.call(cmd_total, shell=True)
    return status


def load_mapped_reads(settings, sample):
    file_in = open(mapped_reads_filename(settings, sample), 'rU')
    file_data = file_in.readlines()
    if len(file_data) > 0:
        return int(file_data[0])
    else:
        return 0


def load_gene_lengths(settings, sample):
    gene_lengths = {}
    data_reader = csv.reader(open(gene_length_filename(settings, sample), 'rU'), delimiter='\t')
    header = next(data_reader)
    for row in data_reader:
        if len(row) == 2:
            gene_lengths[row[0]] = int(row[1])
    return gene_lengths


def read_count_file(filename):
    """ Read the count file generated by HTSeq
        count_data is a dict (tag -> count)
    """
    count_data = {}
    data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
    for row in data_reader:
        # Check that data exists and is not reporting by HTSeq
        if len(row) == 2:
            if len(row[0]) > 2 and row[0][0:2] != '__':
                count_data[row[0]] = int(row[1])
    return count_data


def combine_counts(counts, sample_names):
    """ Combine a set of count dictionaries
        counts is a dictorinary of count_data where key is sample name
    """
    full_tag_list = []
    num_of_samples = len(sample_names)
    # Generate the complete tag list (some samples won't have some tags)
    for sample in sample_names:
        for featuredict in counts[sample]:
            full_tag_list = full_tag_list + [x for x in featuredict]
    full_tag_list = list(set(full_tag_list))
    # Generate matrix zero matrix
    count_matrix = {}
    for tag in full_tag_list:
        count_matrix[tag] = [0] * num_of_samples
    # Update where count exists
    for sample_idx in range(num_of_samples):
        sample = sample_names[sample_idx]
        for featuredict in counts[sample]:
            for tag in featuredict:
                count_matrix[tag][sample_idx] = featuredict[tag] 
    return count_matrix


def save_count_matrix(count_matrix, sample_names, filename):
    """ Save a count_matrix with the sample_names to file
    """
    f_out = open(filename, 'w')
    # Write the header
    f_out.write('\t'.join(['gene_name'] + sample_names) + '\n')
    for tag in sorted(count_matrix):
        count_strs = [str(x) for x in count_matrix[tag]]
        f_out.write('\t'.join([tag] + count_strs) + '\n')
    f_out.close()


def save_mapped_reads_matrix(mapped_reads, sample_names, filename):
    f_out = open(filename, 'w')
    f_out.write('sample\ttotal_mapped_reads\n')
    for s in sample_names:
        f_out.write(s + '\t' + str(mapped_reads[s]) + '\n')
    f_out.close()


def save_gene_length_matrix(gene_lengths, filename):
    f_out = open(filename, 'w')
    f_out.write('gene\tlength\n')
    seen = []
    for s in gene_lengths.keys():
        for gene in gene_lengths[s].keys():
            if gene not in seen:
                f_out.write(gene + '\t' + str(gene_lengths[s][gene]) + '\n')
                seen.append(gene)
    f_out.close()


def count_matrix(settings):
    counts = {}
    for sample in settings.keys():
        if sample != 'None':
            counts[sample] = read_count_file(count_filename(settings, sample))
    sample_names = counts.keys()
    count_matrix = combine_counts(counts, sample_names)
    save_count_matrix(count_matrix, sample_names,
                      settings['None']['output_path'] + 'read_count.matrix')


def gene_lengths(settings, sample, features='gene', attribute='Name'):
    """ Calculate the gene lengths from set of GTF references
    """
    len_file = gene_length_filename(settings, sample)
    f_out = open(len_file, 'w')
    f_out.write('gene_name\tlength\n')
    seen = []
    data_reader = csv.reader(open(settings[sample]['gff_file'], 'rU'), delimiter='\t')
    for row in data_reader:
        for feature in features:
            if len(row) == 9 and row[2] == feature:
                attribs = row[8].split(';')
                for el in attribs:
                    key = el.split('=')[0]
                    value = el.split('=')[1]
                    if key == attribute and value not in seen:
                        gene_length = int(row[4]) - int(row[3]) + 1
                        f_out.write(value + '\t' + str(gene_length) + '\n')
                        seen.append(value)
    f_out.close()
    if os.stat(len_file).st_size == 0:
        return 1
    return 0


def make_profile(settings, sample):
    """ Calculate transcription profile for given regions in BED file
        http://seqanswers.com/forums/showthread.php?t=29399
    """
    if settings[sample]['R2_fastq_file'] == '':
        # https://www.biostars.org/p/14378/
        fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
        cmd_fwd_coverage = 'samtools view -b -F 20 ' + bam_filename(settings, sample,
                                                                    extension=True) + ' > ' + fwd_filename + ' && ' + \
                           'samtools index ' + fwd_filename + ' && ' + \
                           'bedtools coverage -d -abam ' + fwd_filename + ' -b ' + settings[sample]['bed_file'] + \
                           ' > ' + profile_fwd_filename(settings, sample)
        print("Making forward profile: " + cmd_fwd_coverage)
        subprocess.call(cmd_fwd_coverage, shell=True)
        rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
        cmd_rev_coverage = 'samtools view -b -f 16 ' + bam_filename(settings, sample,
                                                                    extension=True) + ' > ' + rev_filename + ' && ' + \
                           'samtools index ' + rev_filename + ' && ' + \
                           'bedtools coverage -d -abam ' + rev_filename + ' -b ' + settings[sample]['bed_file'] + \
                           ' > ' + profile_rev_filename(settings, sample)
        print("Making reverse profile: " + cmd_rev_coverage)
        status = subprocess.call(cmd_rev_coverage, shell=True)

    else:
        fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
        fwd1_filename = bam_filename(settings, sample, extension=False) + '.fwd.1.bam'
        fwd2_filename = bam_filename(settings, sample, extension=False) + '.fwd.2.bam'
        cmd_fwd_coverage = 'samtools view -b -f 83 ' + bam_filename(settings, sample,
                                                                    extension=True) + ' > ' + fwd1_filename + ' && ' + \
                           'samtools index ' + fwd1_filename + ' && ' + \
                           'samtools view -b -f 163 ' + bam_filename(settings, sample,
                                                                     extension=True) + ' > ' + fwd2_filename + ' && ' + \
                           'samtools index ' + fwd2_filename + ' && ' + \
                           'samtools merge -f ' + fwd_filename + ' ' + fwd1_filename + ' ' + fwd2_filename + ' && ' + \
                           'samtools index ' + fwd_filename + ' && ' + \
                           'bedtools coverage -d -abam ' + fwd_filename + ' -b ' + settings[sample]['bed_file'] + \
                           ' > ' + profile_fwd_filename(settings, sample)
        print("Making forward profile: " + cmd_fwd_coverage)
        subprocess.call(cmd_fwd_coverage, shell=True)
        rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
        rev1_filename = bam_filename(settings, sample, extension=False) + '.rev.1.bam'
        rev2_filename = bam_filename(settings, sample, extension=False) + '.rev.2.bam'
        cmd_rev_coverage = 'samtools view -b -f 99 ' + bam_filename(settings, sample,
                                                                    extension=True) + ' > ' + rev1_filename + ' && ' + \
                           'samtools index ' + rev1_filename + ' && ' + \
                           'samtools view -b -f 147 ' + bam_filename(settings, sample,
                                                                     extension=True) + ' > ' + rev2_filename + ' && ' + \
                           'samtools index ' + rev2_filename + ' && ' + \
                           'samtools merge -f ' + rev_filename + ' ' + rev1_filename + ' ' + rev2_filename + ' && ' + \
                           'samtools index ' + rev_filename + ' && ' + \
                           'bedtools coverage -d -abam ' + rev_filename + ' -b ' + settings[sample]['bed_file'] + \
                           ' > ' + profile_rev_filename(settings, sample)
        print("Making reverse profile: " + cmd_rev_coverage)
        status = subprocess.call(cmd_rev_coverage, shell=True)
    return status


def norm_fpkm(settings, bin_path=''):
    """ Calculate normalised RPKM/FPKM expression levels
    """
    cmd_fpkm = 'Rscript ' + bin_path + 'norm_fpkm.r ' + \
               ' ' + count_matrix_filename(settings) + \
               ' ' + mapped_reads_matrix_filename(settings) + \
               ' ' + gene_length_matrix_filename(settings) + \
               ' ' + settings['None']['output_path']
    print("Generating normalised RPKM/FPKMs: " + cmd_fpkm)
    status = subprocess.call(cmd_fpkm, shell=True)
    return status


def de_analysis(settings, group1, group2, output_prefix, bin_path=''):
    """ Calculate DEG analysis between two groups (column numbers 1-indexed)
    """
    cmd_deg = 'Rscript ' + bin_path + 'de_analysis.r ' + \
              ' ' + count_matrix_filename(settings) + \
              ' ' + group1 + \
              ' ' + group2 + \
              ' ' + mapped_reads_matrix_filename(settings) + \
              ' ' + settings['None']['output_path'] + output_prefix
    print("Generating normalised RPKM/FPKMs: " + cmd_deg)
    status = subprocess.call(cmd_deg, shell=True)
    return status


########## CHARACTERIZATION ##########

def load_gff(settings, sample):
    """Function to load gff file data into a dictionary

    Load the information within the gff file into a dictionary. This will be used to generate a list of all parts
    for generation of the context data structure. The first column contains chromosome or plasmid name to serve as main
    keys in the gff dictionary. Mostly the RNAseq analysis focuses on the host plasmid, but information is obtained for
    all entries in the gff.
    The part names are determined by scanning the 9th column in the file for the identifier 'Name='. Parts are further
    classified into different feature types by values in the third column.

    Args:
        filename (str): File path to gff file with all part information

    Return:
        gff (dict): Dictionary with information for every part within the plasmid and host organism

            {
                chromosome:
                    {
                        part1: [part_type, part_direction, start_bp, end_bp, part_attributes],
                        part2: ...
                    }
            }
    """
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
            if part_name is not None:
                if chromo not in gff.keys():
                    gff[chromo] = {}
                gff[chromo][part_name] = [part_type, part_dir, start_bp, end_bp, part_attribs]
    return gff


def load_profiles(settings, sample, normed=False):
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
                new_profile = [cur_start_bp, cur_end_bp, np.zeros(cur_end_bp - cur_start_bp),
                               np.zeros(cur_end_bp - cur_start_bp)]
                new_profile[2][int(row[3]) - 1] = float(row[4])
                profiles[cur_chrom].append(new_profile)
            else:
                cur_profile[0][int(row[3]) - 1] = float(row[4])
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
                cur_profile[1][int(row[3]) - 1] = float(row[4])
    return profiles


def find_profile(profiles, chrom, start_bp, end_bp):
    if chrom in profiles.keys():
        for el in profiles[chrom]:
            if el[0] == start_bp and el[1] == end_bp:
                return [el[2], el[3]]
    return None


def extract_profile_region(profiles, chrom, start_bp, end_bp):
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
                    ext_end_fwd = fwd_profile[:(end_bp - profile_len)]
                    ext_end_rev = rev_profile[:(end_bp - profile_len)]
                new_start_bp = start_bp
                new_end_bp = end_bp
                if ext_start_fwd != []:
                    new_start_bp = 0
                    new_end_bp = end_bp + len(ext_start_fwd)
                new_fwd_profile = ext_start_fwd + fwd_profile + ext_end_fwd
                new_rev_profile = ext_start_rev + rev_profile + ext_end_rev
                region = [new_fwd_profile[new_start_bp:new_end_bp],
                          new_rev_profile[new_start_bp:new_end_bp]]
                break
    return region


def reverse_region(region):
    return [region[1][::-1], region[0][::-1]]


def characterize_promoter_units(settings, sample, upstream_bp=25, downstream_skip_bp=0, downstream_bp=50, normed=False):
    profiles = load_profiles(settings, sample, normed=normed)
    char_data = []
    raw_region = None
    gff = load_gff(settings, sample)
    for chrom in gff.keys():
        for part_name in gff[chrom].keys():
            part_data = gff[chrom][part_name]
            if part_data[0] == 'promoter_unit':
                if part_data[1] == '+':
                    raw_region = extract_profile_region(profiles, chrom,
                                                        (part_data[2] - 1) - upstream_bp,
                                                        part_data[3] + downstream_skip_bp + downstream_bp)
                else:
                    raw_region = reverse_region(extract_profile_region(profiles, chrom,
                                                                       (part_data[
                                                                            2] - 1) - downstream_skip_bp - downstream_bp,
                                                                       part_data[3] + upstream_bp))
                # Calculate performance
                avg_us = np.median(raw_region[0][0:upstream_bp])
                avg_ds = np.median(raw_region[0][-downstream_bp:])
                perf = avg_ds - avg_us
                char_data.append([chrom, part_name, avg_us, avg_ds, perf])
    return char_data


def characterize_promoters(settings, sample, upstream_bp=25, downstream_skip_bp=0, downstream_bp=50, normed=False):
    profiles = load_profiles(settings, sample, normed=normed)
    char_data = []
    raw_region = None
    gff = load_gff(settings, sample)
    for chrom in gff.keys():
        for part_name in gff[chrom].keys():
            part_data = gff[chrom][part_name]
            if part_data[0] == 'promoter':
                if part_data[1] == '+':
                    raw_region = extract_profile_region(profiles, chrom,
                                                        (part_data[2] - 1) - upstream_bp,
                                                        part_data[3] + downstream_skip_bp + downstream_bp)
                else:
                    raw_region = reverse_region(extract_profile_region(profiles, chrom,
                                                                       (part_data[
                                                                            2] - 1) - downstream_skip_bp - downstream_bp,
                                                                       part_data[3] + upstream_bp))
                # Calculate performance
                avg_us = np.median(raw_region[0][0:upstream_bp])
                avg_ds = np.median(raw_region[0][-downstream_bp:])
                perf = avg_ds - avg_us
                char_data.append([chrom, part_name, avg_us, avg_ds, perf])
    return char_data


def characterize_terminators(settings, sample, upstream_bp=50, upstream_skip_bp=0, downstream_bp=25, normed=False):
    profiles = load_profiles(settings, sample, normed=normed)
    char_data = []
    raw_region = None
    gff = load_gff(settings, sample)
    for chrom in gff.keys():
        for part_name in gff[chrom].keys():
            part_data = gff[chrom][part_name]
            if part_data[0] == 'terminator':
                if part_data[1] == '+':
                    raw_region = extract_profile_region(profiles, chrom,
                                                        (part_data[2] - 1) - upstream_skip_bp - upstream_bp,
                                                        part_data[3] + downstream_bp)
                else:
                    raw_region = reverse_region(extract_profile_region(profiles, chrom,
                                                                       (part_data[2] - 1) - downstream_bp,
                                                                       part_data[3] + upstream_skip_bp + upstream_bp))
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
                        t_e = 1.0 / float(avg_us)
                        max_term = 'Y'
                    else:
                        t_e = float(avg_ds) / float(avg_us)
                if t_e != 0.0:
                    t_s = 1.0 / t_e
                else:
                    t_s = -1.0
                char_data.append([chrom, part_name, avg_us, avg_ds, t_e, t_s, max_term])
    return char_data


def characterize_ribozymes(settings, sample, upstream_bp=50, downstream_skip_bp=0, downstream_bp=50, normed=False):
    profiles = load_profiles(settings, sample, normed=normed)
    char_data = []
    raw_region = None
    gff = load_gff(settings, sample)
    for chrom in gff.keys():
        for part_name in gff[chrom].keys():
            part_data = gff[chrom][part_name]
            if part_data[0] == 'ribozyme':
                cut_site = 0
                if 'cut_site' in part_data[4].keys():
                    cut_site = int(part_data[4]['cut_site'])
                if part_data[1] == '+':
                    cur_site_bp = (part_data[2] - 1) + cut_site
                    raw_region = extract_profile_region(profiles, chrom,
                                                        cur_site_bp - upstream_bp,
                                                        cur_site_bp + downstream_skip_bp + downstream_bp)
                else:
                    cur_site_bp = (part_data[3]) - cut_site
                    raw_region = reverse_region(extract_profile_region(profiles, chrom,
                                                                       cur_site_bp - downstream_skip_bp - downstream_bp,
                                                                       cur_site_bp + upstream_bp))
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
                            # c_e = 1.0-(1.0/float(avg_ds))
                            c_e = ((float(avg_ds) + 1.0) / 1.0) - 1.0
                        max_cut = 'Y'
                    else:
                        if avg_ds < avg_us:
                            c_e = 0.0
                        else:
                            # c_e = 1.0-(float(avg_us)/float(avg_ds))
                            c_e = (float(avg_ds) / float(avg_us)) - 1.0
                char_data.append([chrom, part_name, avg_us, avg_ds, c_e, max_cut, cut_site])
    return char_data


def save_characterization_data(settings, sample, data, part_type=None):
    if part_type == 'promoter':
        filename = promoter_profile_perf_filename(settings, sample)
        f_out = open(filename, 'w')
        f_out.write('sample\tchromosome\tpart_name\treads_us\treads_ds\treads_strength\n')
        for d in data:
            f_out.write(sample + '\t' + '\t'.join([str(x) for x in d]) + '\n')
        f_out.close()
    if part_type == 'terminator':
        filename = terminator_profile_perf_filename(settings, sample)
        f_out = open(filename, 'w')
        f_out.write('sample\tchromosome\tpart_name\treads_us\treads_ds\tt_e\tt_s\tmax_term\n')
        for d in data:
            f_out.write(sample + '\t' + '\t'.join([str(x) for x in d]) + '\n')
        f_out.close()
    if part_type == 'ribozyme':
        filename = ribozyme_profile_perf_filename(settings, sample)
        f_out = open(filename, 'w')
        f_out.write('sample\tchromosome\tpart_name\treads_us\treads_ds\tc_e\tmax_cut\tcut_site\n')
        for d in data:
            f_out.write(sample + '\t' + '\t'.join([str(x) for x in d]) + '\n')
        f_out.close()


def combine_promoter_characterizations(settings, samples):
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
            chrom_part_data.append([row[0]] + row[3:])
    f_out = open(combined_promoter_profile_perf_filename(settings), 'w')
    f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\treads_strength\n')
    for chrom in sorted(data.keys()):
        chrom_data = data[chrom]
        for part in sorted(chrom_data.keys()):
            chrom_part_data = chrom_data[part]
            for data_rec in chrom_part_data:
                f_out.write('\t'.join([chrom, part] + data_rec) + '\n')
    f_out.close()
    if os.stat(combined_promoter_profile_perf_filename(settings)).st_size == 0:
        return 1
    return 0


def combine_terminator_characterizations(settings, samples):
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
            chrom_part_data.append([row[0]] + row[3:])
    f_out = open(combined_terminator_profile_perf_filename(settings), 'w')
    f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\tt_e\tt_s\tmax_term\n')
    for chrom in sorted(data.keys()):
        chrom_data = data[chrom]
        for part in sorted(chrom_data.keys()):
            chrom_part_data = chrom_data[part]
            for data_rec in chrom_part_data:
                f_out.write('\t'.join([chrom, part] + data_rec) + '\n')
    f_out.close()
    if os.stat(combined_promoter_profile_perf_filename(settings)).st_size == 0:
        return 1
    return 0


def combine_ribozyme_characterizations(settings, samples):
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
            chrom_part_data.append([row[0]] + row[3:])
    f_out = open(combined_ribozyme_profile_perf_filename(settings), 'w')
    f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\tc_e\tmax_cut\tcut_site\n')
    for chrom in sorted(data.keys()):
        chrom_data = data[chrom]
        for part in sorted(chrom_data.keys()):
            chrom_part_data = chrom_data[part]
            for data_rec in chrom_part_data:
                f_out.write('\t'.join([chrom, part] + data_rec) + '\n')
    f_out.close()
    if os.stat(combined_promoter_profile_perf_filename(settings)).st_size == 0:
        return 1
    return 0


########## FRAGMENT DISTRIBUTIONS ##########

def fragment_length_dists(settings, sample, reads_to_sample=1000000):
    """ Adapted from get_insert_size.py (Wei Li)
    """
    frag_file = fragment_dist_filename(settings, sample)
    sam_file = sam_filename(settings, sample)
    f_in = open(sam_file, 'rU')
    plrdlen = {}
    plrdspan = {}
    objmrl = re.compile('([0-9]+)M$')
    objmtj = re.compile('NH:i:(\d+)')
    nline = 0
    with open(sam_file, 'rU') as ins:
        for lines in ins:
            field = lines.strip().split()
            nline = nline + 1
            if nline >= reads_to_sample:
                break
            if len(field) < 12:
                continue
            try:
                mrl = objmrl.match(field[5])
                if mrl == None:  # ignore non-perfect reads
                    continue
                readlen = int(mrl.group(1))
                if readlen in plrdlen.keys():
                    plrdlen[readlen] = plrdlen[readlen] + 1
                else:
                    plrdlen[readlen] = 1
                if field[6] != '=':
                    continue
                dist = int(field[8])
                if dist <= 0:  # ignore neg dist
                    continue
                mtj = objmtj.search(lines)
                if dist in plrdspan.keys():
                    plrdspan[dist] = plrdspan[dist] + 1
                else:
                    plrdspan[dist] = 1
            except ValueError:
                continue
    f_out = open(frag_file, 'w')
    for k in sorted(plrdspan.keys()):
        f_out.write(str(k) + '\t' + str(plrdspan[k]) + '\n')
    if os.stat(frag_file).st_size == 0:
        return 1
    return 0


def load_norm_factor(settings, sample):
    norm_facs = {}
    norm_fac_file = settings['None']['output_path'] + 'norm.factors.matrix.txt'
    data_reader = csv.reader(open(norm_fac_file, 'rU'), delimiter='\t')
    # Ignore the header
    data_reader.next()
    # Process each line
    for row in data_reader:
        if len(row) == 3:
            norm_facs[row[0]] = [float(row[1]), float(row[2])]
    return (norm_facs[sample][0] * norm_facs[sample][1]) / 1000000.0


def load_fragmentation_dist(settings, sample, max_frag_len=1000):
    frag_dist = np.zeros(max_frag_len + 1)
    frag_file = fragment_dist_filename(settings, sample)
    data_reader = csv.reader(open(frag_file, 'rU'), delimiter='\t')
    # Process each line
    for row in data_reader:
        frag_len = int(row[0])
        frag_count = int(row[1])
        if frag_len <= max_frag_len:
            frag_dist[frag_len] = frag_count
    return frag_dist


def correction_prob(frag_dist, transcript_len, bp_idx):
    frag_sum = 0.0
    for l in range(len(frag_dist)):
        frag_sum += frag_dist[l] * l  # min([l, transcript_len-l+1])
    frag_N = 0.0
    for l in range(len(frag_dist)):
        frag_N += frag_dist[l] * min([l, bp_idx,
                                      (transcript_len - bp_idx + 1),
                                      (transcript_len - l + 1)])
    return frag_N / frag_sum


def correction_profile(frag_dist, transcript_len):
    corr_profile = np.ones(transcript_len)
    for idx in range(transcript_len):
        corr_prob = correction_prob(frag_dist, transcript_len, idx + 1)
        if corr_prob > 0.0:
            corr_profile[idx] = 1.0 / corr_prob
    return corr_profile


def normalise_profiles(settings, sample, correction=False, baseline_us_bp=5):
    """ Will normalise the profiles using TMM factors and correct for edge
        effects using fragment distribution. Some samples have larger RNA
        output, which is unknown in analysis, resulting in genes will be
        undersampled.
    """
    # Load the profiles for the sample
    profiles = load_profiles(settings, sample)
    # If correction of edge effect required
    if correction == True:
        frag_dist = load_fragmentation_dist(settings, sample)
        gff = load_gff(settings, sample)
        for chrom in gff.keys():
            for part_name in gff[chrom].keys():
                part_data = gff[chrom][part_name]
                if part_data[0] == 'transcript':
                    trans_start = part_data[2] - 1
                    trans_end = part_data[3]
                    # Update each profile that contains the transcript
                    for p_data_idx in range(len(profiles[chrom])):
                        p_data = profiles[chrom][p_data_idx]
                        if p_data[0] <= trans_start and p_data[1] >= trans_end:
                            trans_len = trans_end - trans_start

                            # Calc corr_prob_dist (symetrical)
                            corr_prob_dist = correction_profile(frag_dist, trans_len)

                            if part_data[1] == '+':
                                rel_start = trans_start - p_data[0]
                                rel_end = trans_end - p_data[0]
                                # Calc baseline to subtract
                                bl_start = rel_start - baseline_us_bp
                                baseline_reads = 0.0
                                if bl_start >= 0:
                                    baseline_reads = np.mean(p_data[2][bl_start:rel_start])
                                # Calculate profile corrected for edge effects and update
                                corr_profile = corr_prob_dist * (p_data[3][rel_start:rel_end] - baseline_reads).clip(
                                    min=0.0)
                                profiles[chrom][p_data_idx][2][rel_start:rel_end] = corr_profile + baseline_reads
                            else:
                                rel_start = trans_start - p_data[0]
                                rel_end = trans_end - p_data[0]
                                # Calc baseline to subtract
                                bl_start = rel_end + baseline_us_bp
                                baseline_reads = 0.0
                                if bl_start >= p_data[1] - p_data[0]:
                                    baseline_reads = np.mean(p_data[3][rel_end:bl_start])
                                # Calculate profile corrected for edge effects and update
                                corr_profile = corr_prob_dist * (p_data[3][rel_start:rel_end] - baseline_reads).clip(
                                    min=0.0)
                                profiles[chrom][p_data_idx][3][rel_start:rel_end] = corr_profile + baseline_reads

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
                cur_fwd_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx + 1),
                                str(cur_fwd_profile[idx] / norm_fac)]
                cur_rev_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx + 1),
                                str(cur_rev_profile[idx] / norm_fac)]
                f_out_fwd.write('\t'.join(cur_fwd_data) + '\n')
                f_out_rev.write('\t'.join(cur_rev_data) + '\n')
    f_out_fwd.close()
    f_out_rev.close()
    if os.stat(fwd_file).st_size == 0:
        return 1
    if os.stat(rev_file).st_size == 0:
        return 1
    return 0


########## FITTING OF RESPONSE FUNCTIONS ##########

def load_promoter_characterizations(filename, samples):
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


def hill_func(x, Pmin, Pmin_inc, K, n, repress=False):
    if x < 0.0:
        x = 0.0
    if repress == True:
        return Pmin + (Pmin + Pmin_inc) * (math.pow(K, n) / (math.pow(K, n) + math.pow(x, n)))
    else:
        return Pmin + (Pmin + Pmin_inc) * (math.pow(x, n) / (math.pow(K, n) + math.pow(x, n)))


def extract_fit_params(x, P_names, P_types):
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
            fit_params[p_name]['Pmin_inc'] = x[cur_idx + 1]
            cur_idx += 2
        else:
            # If a repessed promoter then 4 parameters
            fit_params[p_name] = {}
            fit_params[p_name]['Pmin'] = x[cur_idx]
            fit_params[p_name]['Pmin_inc'] = x[cur_idx + 1]
            fit_params[p_name]['K'] = x[cur_idx + 2]
            fit_params[p_name]['n'] = x[cur_idx + 3]
            cur_idx += 4
    return fit_params


def promoter_unit_err_func(x, exp_data, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names, P_types,
                           fac_reduce_size):
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
                    fit_outs.append(fit_params[cur_p]['Pmin'] + fit_params[cur_p]['Pmin_inc'])
                else:
                    fit_outs.append(fit_params[cur_p]['Pmin'])
            else:
                # Repressor promoter
                input_val = exp_data[sample][chrom_inputs[p_idx]][PU_inputs[p_idx]]
                fit_outs.append(hill_func(input_val,
                                          fit_params[cur_p]['Pmin'],
                                          fit_params[cur_p]['Pmin_inc'],
                                          fit_params[cur_p]['K'],
                                          fit_params[cur_p]['n'],
                                          repress=True))
        fit_val = np.sum(fit_outs)
        # Numbers are generally in read-depths (reduce to reduce squared values being too large)
        err_diffs.append((exp_val - fit_val) / fac_reduce_size)
    # Return SSE
    return np.sum(np.power(err_diffs, 2.0))


def extract_key_vals(data_str):
    split_str = data_str.split(':')
    key_vals = {}
    for el in split_str:
        key_val_split = el.split('>')
        if len(key_val_split) == 2:
            key_vals[key_val_split[0]] = float(key_val_split[1])
    return key_vals


def extract_list(data_str):
    split_str = data_str.split(',')
    return split_str


def fit_promoter_response_functions(settings, samples, output_name, fac_reduce_size=1000.0):
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
                        x0[cur_param_idx + 1] = 1000.0
                        # Set bounds
                        bnds.append((0.0, None))
                        bnds.append((0.0, None))
                        cur_param_idx += 2
                    else:
                        # Initial conditions (start realistic to improve fitting)
                        x0[cur_param_idx] = 0.0
                        x0[cur_param_idx + 1] = 1000.0
                        x0[cur_param_idx + 2] = 100.0
                        x0[cur_param_idx + 3] = 2.0
                        # Set bounds
                        bnds.append((0.0, None))
                        bnds.append((0.0, None))
                        bnds.append((1.0, None))
                        bnds.append((1.0, 3.9))
                        cur_param_idx += 4
                # methods = BFGS, nelder-mead, Powell, TNC, SLSQP,  L-BFGS-B
                res = scipy.optimize.minimize(promoter_unit_err_func, x0, args=(
                    pro_data, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names, P_types, fac_reduce_size),
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
                        fitted_pro_params[p_name]['Pmin_inc'] = res.x[cur_param_idx + 1]
                        cur_param_idx += 2
                    else:
                        fitted_pro_params[p_name] = {}
                        fitted_pro_params[p_name]['Pmin'] = res.x[cur_param_idx]
                        fitted_pro_params[p_name]['Pmin_inc'] = res.x[cur_param_idx + 1]
                        fitted_pro_params[p_name]['K'] = res.x[cur_param_idx + 2]
                        fitted_pro_params[p_name]['n'] = res.x[cur_param_idx + 3]
                        cur_param_idx += 4
    # Save the results to file
    out_filename = combined_fitted_promoter_perf_filename(settings, output_name)
    f_out = open(out_filename, 'w')
    for p_name in fitted_pro_params.keys():
        f_out.write(p_name)
        for param in fitted_pro_params[p_name].keys():
            f_out.write('\t' + param + '\t' + str(fitted_pro_params[p_name][param]))
        f_out.write('\n')
    f_out.close()


def merge_reus(part_reu_file, profile_perf_filename):
    reu = {}
    with open(part_reu_file) as f:
        datas = f.readlines()
        for line in datas:
            if "terminator" in part_reu_file:
                reu[line.split()[0]] = line.split()[1]
            elif "promoter" in part_reu_file:
                reu[line.split()[0]] = line.split()[3]
        with open(profile_perf_filename.rstrip("txt") + "reu.txt", "w") as out:
            with open(profile_perf_filename) as d:
                datas = csv.reader(d, delimiter='\t')
                for row in datas:
                    if "chromosome" in row:
                        row.append("REU_strength")
                        out.write("\t".join(row) + "\n")
                        continue
                    if "terminator" in profile_perf_filename:
                        part_name = row[1]
                        part_reu = reu[part_name]
                        row.append(part_reu)
                        out.write("\t".join(row) + "\n")
                    elif "promoter" in profile_perf_filename:
                        part_name = row[1].split("-")[1:]
                        if len(part_name) < 2:
                            part = part_name[0]
                            part_reu = reu[part]
                            row.append(part_reu)
                            out.write("\t".join(row) + "\n")
                        else:
                            p1 = part_name[0]
                            p2 = part_name[1]
                            part1_reu = reu[p1]
                            part2_reu = reu[p2]
                            combined_reu = int(part1_reu) + int(part2_reu)
                            reu[row[1]] = str(combined_reu)
                            row.append(str(combined_reu))
                            out.write("\t".join(row) + "\n")
    if os.stat(profile_perf_filename.rstrip("txt") + "reu.txt").st_size == 0:
        print("REU merge failed")
        return 1
    return 0


def comparison_graphs(settings, part_perf_reu_file):
    """
    Terminator strengths from Nature Methods 10, 659 (2013) doi:10.1038/nmeth.2515 (SI pg 31)
    Promoter strengths from pre-publication manuscript TE Gorochowski, YJ Park, C Voigt
    A matrix of graphs will be made for two cases, sample by sample and part by part
    """
    plt.ioff()
    # Plot promoter comparison by sample
    if "promoter" in part_perf_reu_file:
        with open(part_perf_reu_file, "r") as f:
            datas = f.readlines()
            promoters = []
            sample_dic = {}
            for line in datas[1:]:
                sample = line.split()[2]
                seq_strength = line.split()[5]
                reu = line.split()[6]
                promoter = line.split()[1]
                promoters.append(line.split()[1])
                if sample not in sample_dic:
                    sample_dic[sample] = [promoter, [seq_strength, reu]]
                elif promoter not in sample_dic[sample]:
                    sample_dic[sample].append(promoter)
                    sample_dic[sample].append([seq_strength, reu])
            fig1 = plt.figure(figsize=(20, 20))
            # Promoters graphed per sample
            samples = len(sample_dic)
            matrix = math.ceil(math.sqrt(samples))
            subplot_number = 1
            fig1.suptitle("RNASeq Expression Strength vs Part REU Strength Comparisons", fontsize=16, y=1.02)
            for samp in sorted(sample_dic):
                ax1 = fig1.add_subplot(matrix, matrix, subplot_number)
                seq_vals = []
                reu_vals = []
                for i in range(len(sample_dic[samp])):
                    val = sample_dic[samp][i]
                    if type(val) == str:
                        if val in promoters:
                            continue
                    elif type(val) == list:
                        seq_vals.append(val[0])
                        reu_vals.append(val[1])
                y_fpkm = np.array([z for z in seq_vals])
                x_reu = np.array([z for z in reu_vals])
                ax1.set_ylabel("RNA Expression Strength" + "\n" + "(FPKM us / FPKM ds)")
                ax1.set_xlabel("Promoter Strength" + "\n" + "(REU fc)")
                ax1.set_title(samp)
                iteration = enumerate(sample_dic[samp])
                # Need to ensure datapoints are taken from corresponding sample promoter:datapoints
                for idx, obj in iteration:
                    prom = obj
                    XY = next(iteration)[1]
                    ax1.annotate(prom, xy=(XY[1], XY[0]), xytext=(-5, 5), ha='left', textcoords='offset points')
                colorset = cm.rainbow(np.linspace(0, 1, len(x_reu)))
                ax1.scatter(x_reu, y_fpkm, c=colorset)
                subplot_number += 1
                fig1.tight_layout()
        promoter_comparison_filename = settings['None']['output_path'] + "promoter_comparisons_by_sample.png"
        fig1.savefig(promoter_comparison_filename)
    # Plot terminator comparisons by sample
    if "terminator" in part_perf_reu_file:
        with open(part_perf_reu_file, "r") as f:
            datas = f.readlines()
            terminators = []
            sample_dic = {}
            for line in datas[1:]:
                sample = line.split()[2]
                seq_strength = line.split()[6]  # Reads downstream of terminator divided by reads upstream
                reu = line.split()[8]
                terminator = line.split()[1]
                terminators.append(line.split()[1])
                if sample not in sample_dic:
                    sample_dic[sample] = [terminator, [seq_strength, reu]]
                elif terminator not in sample_dic[sample]:
                    sample_dic[sample].append(terminator)
                    sample_dic[sample].append([seq_strength, reu])
            fig3 = plt.figure(figsize=(20, 20))
            samples = len(sample_dic)
            matrix = math.ceil(math.sqrt(samples))
            subplot_number = 1
            fig3.suptitle("RNASeq Expression Strength vs Part REU Strength Comparisons", fontsize=16, y=1.02)
            for samp in sorted(sample_dic):
                ax1 = fig3.add_subplot(matrix, matrix, subplot_number)
                seq_vals = []
                reu_vals = []
                for i in range(len(sample_dic[samp])):
                    val = sample_dic[samp][i]
                    if type(val) == str:
                        if val in terminators:
                            continue
                    elif type(val) == list:
                        seq_vals.append(val[0])
                        reu_vals.append(val[1])
                y_fpkm = np.array([z for z in seq_vals])
                x_reu = np.array([z for z in reu_vals])
                ax1.set_ylabel("RNA Expression Strength" + "\n" + "(FPKM ds / FPKM us)")
                ax1.set_xlabel("Terminator Strength" + "\n" + "(REU fc)")
                ax1.set_title(samp)
                iteration = enumerate(sample_dic[samp])
                # Need to ensure datapoints are taken from corresponding sample promoter:datapoints
                for idx, obj in iteration:
                    term = obj
                    XY = next(iteration)[1]
                    ax1.annotate(term, xy=(XY[1], XY[0]), xytext=(-5, 5), ha='left', textcoords='offset points')
                colorset = cm.rainbow(np.linspace(0, 1, len(x_reu)))
                ax1.scatter(x_reu, y_fpkm, c=colorset)
                subplot_number += 1
                fig3.tight_layout()
        terminator_comparison_filename = settings['None']['output_path'] + "terminator_comparisons_by_sample.png"
        fig3.savefig(terminator_comparison_filename)
    # Plot promoter comparisons by promototer
    if "promoter" in part_perf_reu_file:
        with open(part_perf_reu_file, "r") as f:
            datas = f.readlines()
            promoters = []
            sample_dic = {}
            for line in datas[1:]:
                sample = line.split()[2]
                seq_strength = line.split()[5]
                reu = line.split()[6]
                promoter = line.split()[1]
                promoters.append(line.split()[1])
                if sample not in sample_dic:
                    sample_dic[sample] = [promoter, [seq_strength, reu]]
                elif promoter not in sample_dic[sample]:
                    sample_dic[sample].append(promoter)
                    sample_dic[sample].append([seq_strength, reu])
            fig4 = plt.figure(figsize=(20, 20))
            # Samples graphed per promoter
            proms = len(set(promoters))
            matrix = math.ceil(math.sqrt(proms))
            subplot_number = 1
            fig4.suptitle("RNASeq Expression Strength vs Part REU Strength Comparisons", fontsize=16, y=1.02)
            for prom in sorted(set(promoters)):
                ax1 = fig4.add_subplot(matrix, matrix, subplot_number)
                seq_vals = []
                reu_vals = []
                for i in sample_dic:
                    promoters_and_values = sample_dic[i]  # Only grab promoters and values for this sample
                    idx = promoters_and_values.index(prom)  # Get index for specific promoter within this sample
                    promoter = promoters_and_values[idx]  # Ensure promoter name
                    val = promoters_and_values[idx + 1]  # Ensure values for promoter
                    if type(val) == str:
                        if val in promoters:
                            continue
                    elif type(val) == list:
                        seq_vals.append(val[0])
                        reu_vals.append(val[1])
                y_fpkm = np.array([z for z in seq_vals])
                x_reu = np.array([z for z in reu_vals])
                ax1.set_ylabel("RNA Expression Strength" + "\n" + "(FPKM us / FPKM ds)")
                ax1.set_xlabel("Promoter Strength" + "\n" + "(REU fc)")
                ax1.set_title(prom)
                # Need to ensure datapoints are taken from corresponding sample promoter:datapoints
                k = 0
                for i in sample_dic:
                    XY = [x_reu[k], y_fpkm[k]]
                    ax1.annotate(i, xy=(XY[0], XY[1]), ha='left', xytext=(5, 0), textcoords='offset points')
                    k += 1
                colorset = cm.rainbow(np.linspace(0, 1, len(x_reu)))
                ax1.scatter(x_reu, y_fpkm, c=colorset)
                ax1.ticklabel_format(useOffset=False)
                ax1.set_xlim([int(min(reu_vals)) - 1, int(max(reu_vals)) + 1])
                subplot_number += 1
                fig4.tight_layout()
        promoter_comparison_filename2 = settings['None']['output_path'] + "promoter_comparisons_by_part.png"
        fig4.savefig(promoter_comparison_filename2)
    # PLot terminator comparisons by terminator
    if "terminator" in part_perf_reu_file:
        with open(part_perf_reu_file, "r") as f:
            datas = f.readlines()
            terminators = []
            sample_dic = {}
            for line in datas[1:]:
                sample = line.split()[2]
                seq_strength = line.split()[5]
                reu = line.split()[8]
                terminator = line.split()[1]
                terminators.append(line.split()[1])
                if sample not in sample_dic:
                    sample_dic[sample] = [terminator, [seq_strength, reu]]
                elif terminator not in sample_dic[sample]:
                    sample_dic[sample].append(terminator)
                    sample_dic[sample].append([seq_strength, reu])
            fig5 = plt.figure(figsize=(20, 20))
            # Samples graphed per promoter
            terms = len(set(terminators))
            matrix = math.ceil(math.sqrt(terms))
            subplot_number = 1
            fig5.suptitle("RNASeq Expression Strength vs Part REU Strength Comparisons", fontsize=16, y=1.02)
            for term in sorted(set(terminators)):
                ax1 = fig5.add_subplot(matrix, matrix, subplot_number)
                seq_vals = []
                reu_vals = []
                for i in sample_dic:
                    terminators_and_values = sample_dic[i]  # Only grab promoters and values for this sample
                    idx = terminators_and_values.index(term)  # Get index for specific promoter within this sample
                    terminator = terminators_and_values[idx]  # Ensure promoter name
                    val = terminators_and_values[idx + 1]  # Ensure values for promoter
                    if type(val) == str:
                        if val in terminators:
                            continue
                    elif type(val) == list:
                        seq_vals.append(val[0])
                        reu_vals.append(val[1])
                y_fpkm = np.array([z for z in seq_vals])
                x_reu = np.array([z for z in reu_vals])
                ax1.set_ylabel("RNA Expression Strength" + "\n" + "(FPKM us / FPKM ds)")
                ax1.set_xlabel("Terminator Strength" + "\n" + "(REU fc)")
                ax1.set_title(term)
                # Need to ensure datapoints are taken from corresponding sample promoter:datapoints
                k = 0
                for i in sample_dic:
                    XY = [x_reu[k], y_fpkm[k]]
                    ax1.annotate(i, xy=(XY[0], XY[1]), ha='left', xytext=(5, 0), textcoords='offset points')
                    k += 1
                colorset = cm.rainbow(np.linspace(0, 1, len(x_reu)))
                ax1.scatter(x_reu, y_fpkm, c=colorset)
                ax1.ticklabel_format(useOffset=False)
                ax1.set_xlim([int(min(reu_vals)) - 1, int(max(reu_vals)) + 1])
                subplot_number += 1
                fig5.tight_layout()
        terminator_comparison_filename2 = settings['None']['output_path'] + "terminator_comparisons_by_part.png"
        fig5.savefig(terminator_comparison_filename2)
    return 0


########## VARIANT ANALYSIS ##########

def mark_dupes(settings, sample):
    """Mark duplicates using Picard-tools
    """
    cmd_markdupes = 'java -Xmx8g -jar /seq/software/picard/current/bin/picard.jar MarkDuplicates' + \
      ' I=' + bam_filename(settings, sample, extension=True) + \
      ' O=' + bam_duplicate_name(settings, sample, extension=True) + \
      ' CREATE_INDEX=true' + \
      ' VALIDATION_STRINGENCY=SILENT' + \
      ' M=' + settings[sample]['output_path'] + 'marked_dup_metrics.txt'

    print("Picard mark duplicates: " + cmd_markdupes)
    status = subprocess.call(cmd_markdupes, shell=True)
    if status != 0:
        return 1
    return status


def make_ref_dictionary(settings, sample):
    """Make reference dictionary with Picard-tools
    """
    cmd_refdict = 'java -Xmx1g -jar /seq/software/picard/current/bin/picard.jar CreateSequenceDictionary' + \
    ' R=' + settings[sample]['fasta_file'] + \
    ' O=' + settings[sample]['fasta_file'].rstrip(".fasta") + '.dict'

    print("Making reference dictionary: " + cmd_refdict)
    status = subprocess.call(cmd_refdict, shell=True)
    if status != 0:
        return 1
    return status


def add_readgroups(settings, sample):
    """Add readgroups to marked dupe bam file with Picard-tools
    """
    cmd_addreadgroups = 'java -Xmx4g -jar /seq/software/picard/current/bin/picard.jar AddOrReplaceReadGroups' + \
    ' I=' + bam_duplicate_name(settings, sample, extension=True) + \
    ' O=' + bam_readgroup_name(settings, sample, extension=True) + \
    ' SO=coordinate RGID={id} RGLB=foundry RGPL=illumina RGPU={id} RGSM={sample} && samtools index {rgbam}'.format(
        sample = sample, id = sample.split("_")[-1], rgbam = bam_readgroup_name(settings, sample, extension=True))
    print("Adding readgroups to marked dupe bam: " + cmd_addreadgroups)
    status = subprocess.call(cmd_addreadgroups, shell=True)
    if status != 0:
        return 1
    return status


def gatk_haplotype_caller(settings, sample):
    """Function to call variants with GATK Haplotype Caller (HC)

    GATK HaplotypeCaller is used to collect vcf information from a mapping file
        -L name of the circuit design in the reference fasta
        -R reference fasta file, can contain circuit & host
        -I input bam file that has read groups added, duplicates marked
        -stand_call_conf 20 and emit_conf 20 set the threshold for discarding variants to Q20. Less stringent since it
         is RNAseq
        -o output file name sample_#.vcf

    Args:
        settings: (dict) The settings dictionary as described in load_settings()
        samples: (str) A sample to generate a VCF file for, one sample is operated on per HC run

    Returns:
        status: (int) If subprocess call succeeds return 0, else return 1 for failure
    """
    cmd_haplotyper = 'java -Xmx4g -jar /broad/software/free/Linux/redhat_6_x86_64/pkgs/GATK3-2.2/GenomeAnalysisTK.jar' + \
    ' -T HaplotypeCaller' + \
    ' -R ' + settings[sample]['fasta_file'] + \
    ' -L ' + os.path.basename(settings[sample]['fasta_file']).rstrip('.fasta') + \
    ' -I ' + bam_readgroup_name(settings, sample, extension=True) + \
    ' -stand_call_conf 20.0' + \
    ' -stand_emit_conf 20.0' + \
    ' -o ' + vcf_name(settings, sample)

    print("GATK3 HaplotypeCaller: " + cmd_haplotyper)
    status = subprocess.call(cmd_haplotyper, shell=True)
    if status != 0:
        return 1
    return status


########## CONTEXT ANALYSIS ##########

def determine_context(gff):
    """Function to determine the context of a part, looking at up to 2 parts before/after and constructing a schema for
        interpreting the context at a glance

    The list is sorted by position so the first part will be at the start of the plasmid, and the last part will be at
    the end. This can be traversed for promoters, terminators, ribozymes, and genes while constructing appropriate
    contexts based upon part combinations.
    For a new part type other than promoter, ribozyme, gene, terminator, an explicit new device abreviation needs to
    be created.
    The context for each part is constructed to highlight the part in '' and surrounding parts with -.
        e.g. Promoter, promoter, ribozyme, gene. Ribozyme is the part queried, thus expressed as:
                P-P'r'CDS

    Args:
        gff (dict): gff file of circuit

            {
                chrom:
                    {
                        part1: [start bp, end bp, orientation, {info}
                        part2: ...
                    }
            }

    Returns:
        Dictionary containing each part in the circuit with 3 to 4 part surrounding context. If more than one context is
        found for a part it is added as another entry in the list of contexts for that part.
            {
                part1: [context1, context2, ...]
                part2: [context1, context2, ...],
            }
    """
    partlist = []
    construct = ''
    for chrom in gff:
        for part in gff[chrom]:
            if gff[chrom][part][0] == 'promoter':
                construct = chrom
                break
    for part, data in gff[construct].items():
            partlist.append([part, data])
    partlist = [x for x in partlist if 'promoter_unit' != x[1][0] and 'transcript' != x[1][0]]
    sorted_parts = sorted(partlist, key=lambda j: j[1][3])
    devices = []
    for part in sorted_parts:
        if part[1][0] == 'promoter':
            devices.append('P')
            continue
        if part[1][0] == 'ribozyme':
            devices.append('r')
            continue
        if part[1][0] == 'gene':
            devices.append('CDS')
            continue
        if part[1][0] == 'terminator':
            devices.append('T')
            continue
    design_context = {}
    for idx, entry in enumerate(sorted_parts):
        if idx == 0:
            design_context[entry[0]] = [devices[idx], devices[idx] + "'" + devices[idx+1] + '-' + devices[idx+2] + '-' +
                                        devices[idx+3]]
        elif idx == 1:
            design_context[entry[0]] = [devices[idx], devices[idx-1] + "'" + devices[idx] + "'" + devices[idx+1] + '-' +
                                        devices[idx+2]]
        else:
            design_context[entry[0]] = [devices[idx], devices[idx-2] + '-' + devices[idx-1] + "'" + devices[idx] + "'" +
                                        '-'.join(devices[idx+1:idx+2])]
    return design_context


def join_gene_data(gff, part_fpkm, prom_perf, term_perf, ribo_perf, context, counts, snp_dic):
    """Funciton to join the different context data dictionaries to a central data structure

    Parts need follow strictly the formatting specified when first initiating the GFF file. Any re-used parts will be
    distinguished by a hyphen and an incremental number. All promoter units will be proceeded by 'PU-', trancsription
    units by 'TR-'.

    Args:
        gff (dict): As described in load_gff()
        part_fpkm (dict): As described in load_fpkm()
        prom_perf (dict): As described in load_perf()
        term_perf (dict): As described in load_perf()
        ribo_perf (dict): As described in load_perf()
        context (dict): As described in determine_context()
        counts (dict): Described in load_counts()
        snp_dic (dict): VCF information, SNPs, indels for each sample

            {
                chrom
                    {
                        position of snp
                            {
                                sample1 : [Ref, Alt, Sum of #read quality score, Allele freuency, Depth],
                                sample2 : ...
                            }
                    }
            }

    Returns:
        gene_info (dict): Joined data of the input dictionaries in new structure

            {
                part name
                    {
                        'name'     : part name all obtained from gff input,
                        'context'  : context of part from context input,
                        'snp'      :
                            {
                                number of de_analysis file  :
                                {
                                    loc : snp_dic[loc], populated if all samples share snp position from input snp_dict,
                                }
                            }
                        'perf'     : based on type of part as specified in gff file, will populate performance info from
                                        input dict,
                        'fpkm'     : fpkm performance from input fpkm dict, reported as fold chance with p-value and FDR
                    }
            }

    Raises:
        KeyError: The part did not reach fold change threshold of 1.8 or -1.8 to be considered functionally impacted and
            was excluded in the fpkm_dict passed to this function, so while iterating through keys (parts) from gff file
            there will be missing keys corresponding to < abs(1.8) fold change parts
    """
    promoter_types = ['P-CDS', 'T-P-CDS', 'CDS-P-CDS', 'P-P', 'P-P-P', 'CDS-P-P', 'T-P-P', 'P-P-CDS', 'P-P-r-CDS',
                      'r-P-P-CDS']
    terminator_types = ['T-P', 'CDS-T-P', 'CDS-T-T', 'T-T-P', 'CDS-T-r-T', 'CDS-T-CDS', 'CDS-T-r-CDS', 'T-r-CDS']
    context_types = set(promoter_types + terminator_types)

    gene_info = {}
    for chrom in prom_perf:
        snplocs = [int(x) for x in snp_dic[chrom].keys()]
        reduced_cxt_set = [k for k, v in gff[chrom].items() if v[0] != 'promoter_unit' and v[0] != 'transcript']
        for part, data in gff[chrom].items():
            if part not in gene_info:
                gene_info[part] = OrderedDict([('name', part),
                                              ('context', set()),
                                              ('snp', []),
                                              ('perf', []),
                                              ('fpkm', {}),
                                              ('counts', {})])
            # Assemble context data
            if part in reduced_cxt_set:
                gene_info[part]['context'].add(context[part][1])
            # Assemble fpkm data
            for de in range(len(part_fpkm)):
                if de not in gene_info[part]['fpkm']:
                    gene_info[part]['fpkm'][de] = []
                    try:
                        gene_info[part]['fpkm'][de].append(part_fpkm[de][part])
                    except KeyError:
                        pass
            # Assemble snp data
            for loc in snplocs:
                # Allow 10 bp up and downstream for snp identification
                if (min(data[2:4])-10) <= loc <= (max(data[2:4])+10):
                    if loc not in gene_info[part]['snp']:
                        gene_info[part]['snp'] = {}
                    gene_info[part]['snp'][loc] = snp_dic[chrom][str(loc)]
            # Assemble count data
            gene_info[part]['counts'] = counts[part]
            # Assemble performance data
            if part in prom_perf[chrom].keys():
                gene_info[part]['perf'] = prom_perf[chrom][part]
            if part in term_perf[chrom].keys():
                gene_info[part]['perf'] = term_perf[chrom][part]
            if part in ribo_perf[chrom].keys():
                gene_info[part]['perf'] = ribo_perf[chrom][part]
    return gene_info


def load_fpkm(fpkm_file):
    """Funciton to load fpkm data into dictionary from file

    FPKM data is collected for parts if the fold change is greater than 1.8. If more than just the synthetic design FPKM
    is wanted, remove the condition if statement filtering out 'locusTag' data.

    Args:
        fpkm_file: (str) A file path to the fpkm file that was created from de_analysis step of pipeline. The filename
         for de_analysis if dnyamic, but input here is generalized to sample1_vs_sample2 name schema.

    Returns:
        part_fpkm (dict): a dictionary with FPKM fold change, p-value, and FDR (false discovery rate) as calculated
         by this pipeline

        {
            part1 name: [fold change, p-value, FDR],
            part2 name: [fold change, p-value, FDR],
            etc.
        }

    """
    with open(fpkm_file) as f:
        data = f.readlines()
    part_fpkm = {}
    for line in data[1:]:
        cols = line.split()
        if 'locusTag' in cols[0]:
            continue
        if float(cols[1]) >= 1.8 or float(cols[1]) <= -1.8:
            part_fpkm[cols[0]] = ['%.3f' % float(cols[1]), '%.3f' % float(cols[3]), '%.3f' % float(cols[4])]
    return part_fpkm


def load_vcf(vcf_file, sample):
    """Function to load vcf data from one sample file into dictionary

    Args:
        vcf_file: (str) Relative path to vcf files produced from RNAseq pipeline
        sample: (str) Name of sample for which vcf data file will be parsed

    Returns:
        snps: (dict) Data structure containing all snps found in vcf file for a particular sample

            {
                sample: [ [line from vcf file with SNP1 data],
                          [line from vcf file with SNP2 data],
                          etc,
                        ]
            }
    """
    with open(vcf_file) as f:
        data = f.readlines()
    snps = {sample: []}
    for line in data:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            continue
        else:
            cols = line.split()
            snps[sample].append(cols)
    return snps


def load_perf(perf_file):
    """Function to load part performance data, able to load the three types of performances calculated in pipeline

    The parts promoter, terminator, and ribozyme have central control functions over the expression of the genes that
    control the function of either circuits or gene clusters. So these performances are modeled in the genetic_analyzer
    script to assess their performance. Their performance is collected here by parsing those performance files and
    generating a dictionary

    Args:
        perf_file: (str) Relative path to promoter, terminator, or ribozyme performance data

    Returns:
        perfs: (dict) Data structure containing all the parts within each perf file. These will correlate with parts
         defined in the gff file

            {
                chrom :
                    {
                        part :
                            {
                                sample: [part_strength],
                                sample2: [part_strength]
                                etc.
                            }
                    }
            }
    """
    perfs = {}
    with open(perf_file) as f:
        data = f.readlines()
    for line in data[1:]:
        line = line.split()
        chrom = line[0]
        pu = line[1]
        sample = line[2]
        if chrom not in perfs:
            perfs[chrom] = {}
        if pu not in perfs[chrom]:
            perfs[chrom][pu] = {}
        if sample not in perfs[chrom][pu]:
            if 'promoter' in perf_file:
                perfs[chrom][pu][sample] = line[5]
            elif 'terminator' in perf_file:
                perfs[chrom][pu][sample] = line[6]
            elif 'ribozyme' in perf_file:
                perfs[chrom][pu][sample] = line[5]
    return perfs


def load_counts(count_file):
    """Function to convert fpkm.normed.matrix.txt to a dictionary

    Args:
        count_file: (string) File generated by edgeR normaliziation methods containing FPKM of every gene per sample

    Returns:
        counts: (dict) Data structure that is organized to be accessed by the context effect data join routine.

            {
                part (str) :
                    {
                        sample (str) : (str(float(3 digit precision))
                    }
            }
    """
    with open(count_file) as f:
        datas = f.readlines()
    samples = datas[0].split()
    counts = {}
    for line in datas[1:]:
        part = line.split()[0]
        fpkms = line.split()[1:]
        if part not in counts:
            counts[part] = {}
        for sample in samples:
            if sample not in counts[part]:
                counts[part][sample] = []
            idx = samples.index(sample)
            counts[part][sample] = '%.3f' % float(fpkms[idx])
    return counts


def combine_vcf_data(settings, samples, gff):
    """Function to collect snp data from a VCF file for one sample

    A VCF file is parsed for all snps reported, no filtering is used at this stage. The GFF file is parsed for which
    chromosomes are present in the experiment, however only the host plasmid gets VCF calling performed, so the dict
    will contain values for just the host plasmid chromosome.
    This will not be the case if a synthetic gene cluster is incorporated into a host chromosome. This case will have
    the host chromosome SNP data analyzed prior and contained in this data struct as well.

    Args:
        samples: (str) A sample name corresponding to samples in settings dict, used to key snp dict
        gff: (dict) GFF dictionary containing the name of chromosomes.

    Returns:
        snp_dic: (dict) Data structure for each sample for each SNP position

            {
                chromosome :
                    {
                        snp position :
                            {
                                sample : [Ref allele, alt allele, Sum of QScore for depth, Allele frequency, Depth],
                            }
                    }
            }
    """
    snp_dic = {}
    for sample in samples:
        if sample != 'None':
            try:
                sample_vcf = load_vcf(vcf_name(settings, sample), sample)
            except OSError:
                continue
            for chrom in gff:
                if chrom not in snp_dic:
                    snp_dic[chrom] = {}
                for per_sample_snp in sample_vcf[sample]:
                    cols = [per_sample_snp][0]
                    snp = cols[1]
                    AF = [x for x in cols[7].split(';') if 'AF' == x.split('=')[0]][0]
                    AF = AF[:7]
                    DP = [y for y in cols[7].split(';') if 'DP' == y.split('=')[0]][0]
                    if snp not in snp_dic[chrom]:
                        snp_dic[chrom][snp] = {}
                    if sample not in snp_dic[chrom][snp]:
                        snp_dic[chrom][snp][sample] = []
                    snp_dic[chrom][snp][sample].append(cols[3:6])
                    snp_dic[chrom][snp][sample][0].extend((AF, DP))
    return snp_dic


def save_context_data(settings, sample_part_data):
    """Funciton to save the context effect data to file

    File name is determined from the settings output path for 'None' sample, the folder above sample-level folders.
    Data output format is encapsulated here, any change to format here will be reflected in context data output file.

    Args:
        settings: (dict) As described in load_settings() function
        sample_part_data: The context data dictionary formatted in context_effects.py()

    Returns:
        (int): Status of operation, 0 if successfully creates file after checking for file size. Return 1 for failure
         if file contains no data
    """
    cxt_file = context_filename(settings)
    with open(cxt_file, 'w') as f:
        for x in sample_part_data['feature']:
            f.write('FEATURE' + '\t' + x + '\n')
            for k in sample_part_data['feature'][x]:
                for g, h in k.items():
                    if not h:
                        continue
                    else:
                        if g == 'name':
                            f.write('\t' + g + '\t' + k[g] + '\n')
                        elif g == 'context':
                            f.write('\t\t' + g + '\t\t')
                            f.write([str(e) for e in k[g]][0] + '\n')
                        elif g == 'perf':
                            f.write('\t\t' + 'performance' + '\n')
                            for j, v in k[g].items():
                                f.write('\t\t\t' + j + ": " + v + '\n')
                        elif g == 'fpkm':
                            for de in k[g]:
                                if len(k[g][de]) > 0:
                                    f.write('\t\tgroup1_vs_group2_' + str(de) + '\n')
                                    fpkm, pval, fdr = k[g][de][0]
                                    f.write('\t\t\tfoldchange: ' + str(fpkm) + '\tp-val: ' +pval + '\tfdr: ' + fdr + '\n')
                        elif g == 'snp':
                            f.write('\t\t' + g + '\n')
                            for snp in k[g]:
                                f.write('\t\t\t' + str(snp) + '\t\tsample\t\tAlt\t\tRef\t\tQual\t\tAlleleFreq\tDepth\t\tCounts\n')
                                samplesnps = k[g][snp]
                                for sample, data in samplesnps.items():
                                    f.write('\t\t\t\t\t' + sample + '\t\t' + '\t\t'.join(data[0]) + '\t\t' + k['counts'][sample] + '\n')
        f.close()
    if os.stat(cxt_file).st_size == 0:
        return 1
    return 0
