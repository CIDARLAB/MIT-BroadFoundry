"""
Provide:
Path to directory with GFF, BED, Fasta directorie and files
Path to directory with fastq files
Sample name
"""

import os
import sys
import argparse


def parseargs():
    parser = argparse.ArgumentParser(description="Create settings file for RNASeq of Genes and Parts Pipeline")
    arggroup = parser.add_argument_group("Arguments")
    arggroup.add_argument('--datadir', '-d', action="store", default=None, dest="data_dir",
                          help="Directory with bed, fasta, and gff files " + \
                                "eg. /btl/foundry/data/reference/")
    arggroup.add_argument('--fastqdir', '-f', action="store", default=None, dest="fastq_dir",
                          help="Directory with Fastq files")
    arggroup.add_argument('--names', '-n', action="append", default=None, dest="samplenames",
                          help="Prefix for sample output, 1 or 2 names")
    arggroup.add_argument('--reference', '-r', action="store", default=None, dest="reference",
                          help="Name of reference in /btl/foundry/data/reference/" + \
                                "eg. 0x58v50")
    options = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return options


def main(options):
    data_dir = options.data_dir
    fastq_dir = options.fastq_dir
    samplenames = options.samplenames
    reference = options.reference
    fastq_list = []
    for f in os.listdir(fastq_dir):
        if f.endswith("fastq"):
            fastq_list.append(f)

    fasta = data_dir + "fasta/" + reference + ".fasta"
    gff = data_dir + "gff/" + reference + ".gff"
    bed = data_dir + "bed/" + reference + ".bed"
    if len(samplenames) > 1:
        temp_path = ["./alignment/" + samplenames[0] + "_", "./alignment/" + samplenames[1] + "_"]
        output_path = ["./results/" + samplenames[0] + "_", "./results/" + samplenames[1] + "_"]
    else:
        temp_path = "./alignment/" + samplenames[0] + "_"
        output_path = "./results/" + samplenames[0] + "_"

    header = ["sample", "fasta_file", "gff_file", "bed_file", "R1_fastq_file", "R2_fastq_file", "temp_path",
              "output_path"]

    # if os.path.exists("settings.txt"):
    # sys.exit("Settings file already exists! Move or change the filename before running this script")

    # H9FULADXX.ATTATGTT.unmapped.2_Circuit_sample_T4_2.fastq
    forwards = []
    reverse = []
    for x in fastq_list:
        end = os.path.basename(x).split("_")[0]
        end = end.split(".")[-1]
        if end == '1':
            forwards.append(x)
        elif end == '2':
            reverse.append(x)

    forwards = sorted(forwards)
    reverse = sorted(reverse)
    fastq_dirname = os.path.dirname(fastq_dir)
    with open("settings.txt", 'w') as out:
        out.write("\t".join(x for x in header) + "\n")
        out.write("None\t\t\t\t\t\t\t./results/\n")
        if len(samplenames) == 2:
            # Samples = 1/2 fastq_num
            fastq_num = 0
            for sample_num in range(len(fastq_list) / 4):
                output = [samplenames[0] + "_" + str(sample_num), fasta, gff, bed,
                          fastq_dirname + "/" + forwards[fastq_num], fastq_dirname + "/" + reverse[fastq_num],
                          temp_path[0] + str(sample_num) + "/", output_path[0] + str(sample_num) + "/"]
                out.write("\t".join([x for x in output]))
                out.write("\n")
                output = [samplenames[1] + "_" + str(sample_num), fasta, gff, bed,
                          fastq_dirname + "/" + forwards[fastq_num + 1], fastq_dirname + "/" + reverse[fastq_num + 1],
                          temp_path[1] + str(sample_num) + "/", output_path[1] + str(sample_num) + "/"]
                out.write("\t".join([x for x in output]))
                out.write("\n")
                fastq_num += 2
        else:
            for sample_num in range(len(fastq_list) / 2):
                output = [samplenames[0] + "_" + str(sample_num), fasta, gff, bed,
                          fastq_dirname + "/" + forwards[fastq_num], fastq_dirname + "/" + reverse[fastq_num],
                          temp_path + str(sample_num) + "/", output_path + str(sample_num) + "/"]
                out.write("\t".join([x for x in output]))
                out.write("\n")


if __name__ == "__main__":
    args = parseargs()
    status = main(args)
    exit(status)
