
# From the regular expression of the entire plasmid containing Ns, make 200 bp regular expressions for each read
python ../../bin/regex_from_seq.py ./data/pcr_seqs/HlyIIR_NOR2_XZinv.fa 200 ./data/regex_refs/HlyIIR_2NOR_regexs.txt

# Run the dialout script, see script for arguments description
bsub -e /btl/foundry/users/tom/dialout/runs/2015-12_raissa/results/HlyIIR_2NOR_stderr.txt -o /btl/foundry/users/tom/dialout/runs/2015-12_raissa/results/HlyIIR_2NOR_stdout.txt -q forest -N  "python /btl/foundry/users/tom/dialout/bin/perfect_dialout.py /btl/foundry/users/tom/dialout/runs/2015-12_raissa/data/regex_refs/HlyIIR_2NOR_regexs.txt /btl/projects/Foundry/Lauren/barcodes_Sept2015/spl1_S1_L001_R1_001_test.fastq /btl/projects/Foundry/Lauren/barcodes_Sept2015/spl1_S1_L001_R2_001_test.fastq 20 20 2 3 1 /btl/foundry/users/tom/dialout/runs/2015-12_raissa/results/HlyIIR_2NOR_"
