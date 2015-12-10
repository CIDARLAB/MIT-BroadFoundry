
# Example of the dialout scripts

mkdir results

# From the regular expression of the entire plasmid containing Ns (found in /data/pcr_seqs), this script takes 200 bp from each end of the sequence and generates regular expressions for Ns. These are the pair of regular expressions used to check the design.
python ../../bin/regex_from_seq.py ./data/pcr_seqs/HlyIIR_NOR2_XZinv.fa 200 ./data/regex_refs/HlyIIR_2NOR_regexs.txt

# Run the dialout script, see script for arguments description:
# 1. Location of the regular expression file (for each design 2 regular expressions for each read)
# 2. FASTQ read 1 file
# 3. FASTQ read 2 file
# 4. Length of shared primer for read 1 for all designs (used to check orientation of read)
# 5. Length of shared primer for read 2 for all designs (used to check orientation of read)
# 6. Location (index) of the forward barcode used for extraction (indexes start at 1 for read 1 and -1 for read 2)
# 7. Location (index) of the reverse barcode used for extraction (indexes start at 1 for read 1 and -1 for read 2)
# 8. Other barcode indexes to output (but not used in checking for uniqueness)
# 9. Output path
bsub -e /btl/foundry/users/tom/pool_dialout/runs/2015-12_raissa/results/HlyIIR_2NOR_stderr.txt -o /btl/foundry/users/tom/pool_dialout/runs/2015-12_raissa/results/HlyIIR_2NOR_stdout.txt -q forest -N  "python /btl/foundry/users/tom/pool_dialout/bin/perfect_dialout.py /btl/foundry/users/tom/pool_dialout/runs/2015-12_raissa/data/regex_refs/HlyIIR_2NOR_regexs.txt /btl/foundry/users/tom/pool_dialout/data/spl1_S1_L001_R1_001_sample.fastq /btl/foundry/users/tom/pool_dialout/data/spl1_S1_L001_R2_001_sample.fastq 20 20 2 3 1 /btl/foundry/users/tom/pool_dialout/runs/2015-12_raissa/results/HlyIIR_2NOR_"
