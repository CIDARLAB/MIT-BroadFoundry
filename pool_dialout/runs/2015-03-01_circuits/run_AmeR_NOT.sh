qsub -N AmeR_NOT -e /home/unix/tgorocho/pool_mapping/runs/2015-03-01_circuits/results/AmeR_NOT_stderr.txt -o /home/unix/tgorocho/pool_mapping/runs/2015-03-01_circuits/results/AmeR_NOT_stdout.txt -b y -q gaag "python /home/unix/tgorocho/pool_mapping/bin/perfect_dialout.py /home/unix/tgorocho/pool_mapping/runs/2015-03-01_circuits/data/regex_refs/AmeR_NOT_regexs.txt /btl/data/MiSeq0/runs/lauren/150227_M03102_0065_000000000-AE6TR/Data/Intensities/BaseCalls/AmeR-NOT-tagged-57_S1_L001_R1_001.fastq /btl/data/MiSeq0/runs/lauren/150227_M03102_0065_000000000-AE6TR/Data/Intensities/BaseCalls/AmeR-NOT-tagged-57_S1_L001_R2_001.fastq 20 20 2 3 /home/unix/tgorocho/pool_mapping/runs/2015-03-01_circuits/results/AmeR_NOT_"