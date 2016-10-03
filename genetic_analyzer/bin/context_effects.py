#!/usr/bin/env python2.7

#   Alex Cristofaro Broad Foundry
#   OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

"""Context effect collector from RNAseq pipeline data

This module will populate context effects for all parts described in the input GFF file to RNAseq pipeline.
Two major output types are grouped into Part and Feature effects. The part data contains all context effects on a per-
part basis. The feature data contains all context effects grouped into the part types, such as promoters, terminators,
ribozymes, genes, promoter_units, and transcript.

Example:
    $ python context_effects.py -settings settings.txt

Todo:
    * Include state from GFF file and inferencing of success status per sample
    * Handle more than one DGE scenario
"""
import sys, os
import genetic_analyzer as ga
import argparse


def main():
    """Main function of module, will load all required data structures from various files to construct a central
        data structure containing all context effect data

    FPKMs are obtained from edgeR TMM DE gene calculations after obtaining counts from HTSeq
    Context is derived from the GFF file.
    Snps come from GATK haplotype caller VCF.
    Performance profiles are calculated from genetic_analyzer functions after looking at bedtools coverage data
    Normalized FPKM counts are displayed alongside any SNP data from fpkm.normed.matrix.txt

    Returns:
        status: (int) If file is found to contain data after writing finishes, success exits with status 0. Else
                failure status is indicated by exit status 1
    """
    # Parse the command line inputs
    parser = argparse.ArgumentParser(description="context effects")
    parser.add_argument("-settings", dest="settings", required=True, help="settings.txt", metavar="string")
    args = parser.parse_args()
    settings = ga.load_settings(args.settings)
    samples = [s for s in settings if s != 'None']
    resultsdir = os.listdir(settings['None']['output_path'])
    gff = ga.load_gff(settings, samples[0])
    fpkm_files = [os.path.abspath(settings['None']['output_path'] + x) for x in resultsdir if 'de.analysis.txt' in x]
    part_fpkm = []
    for f in fpkm_files:
        part_fpkm.append(ga.load_fpkm(f))
    prom_perf = ga.load_perf(ga.combined_promoter_profile_perf_filename(settings))
    term_perf = ga.load_perf(ga.combined_terminator_profile_perf_filename(settings))
    ribo_perf = ga.load_perf(ga.combined_ribozyme_profile_perf_filename(settings))

    context = ga.determine_context(gff)

    snp_dic = ga.combine_vcf_data(settings, samples, gff)

    counts = ga.load_counts(ga.normed_counts_filename(settings))
    
    gene_info = ga.join_gene_data(gff, part_fpkm, prom_perf, term_perf, ribo_perf, context, counts, snp_dic)

    for sample in settings:
        if sample != 'None':
            sample_part_data = {'part': {},
                                'feature': {
                                       'transcript': [],
                                       'promoter': [],
                                       'promoter_unit': [],
                                       'terminator': [],
                                       'ribozyme': [],
                                       'gene': []
                                            }
                                }
            for chrom in prom_perf:
                # Todo: Need better way of reducing analysis set to host plasmid instead of every \
                # Todo: chromosome in GFF file. Querying only the chromosome in perf is a 'hack'
                for part in gff[chrom]:
                    sample_part_data['part'][part] = gene_info[part]
                    for feature in sample_part_data['feature']:
                        if gff[chrom][part][0] == feature:
                            if part not in sample_part_data['feature'][feature]:
                                try:
                                    sample_part_data['feature'][feature].append(gene_info[part])
                                except KeyError:
                                    pass

    status = ga.save_context_data(settings, sample_part_data)
    return status


if __name__ == "__main__":
    status_main = main()
    sys.exit(status_main)
