bsub -e /home/unix/tgorocho/pool_mapping/runs/2015-09_nif/results/nifK_stderr.txt -o /home/unix/tgorocho/pool_mapping/runs/2015-09_nif/results/nifK_stdout.txt -q forest -N  "python /home/unix/tgorocho/pool_mapping/bin/perfect_dialout.py /home/unix/tgorocho/pool_mapping/runs/2015-09_nif/data/regex_refs/nifK_regexs.txt /btl/projects/Foundry/Yongjin/Tom/nif_fastq/nifK_R1.fastq /btl/projects/Foundry/Yongjin/Tom/nif_fastq/nifK_R2.fastq 20 20 1 2 /home/unix/tgorocho/pool_mapping/runs/2015-09_nif/results/nifK_"