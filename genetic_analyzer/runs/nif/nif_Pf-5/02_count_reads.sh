
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Count mapped reads in gene features

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/02_count_reads_N3.out.log -e ./logs/02_count_reads_N3.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N3 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N4.out.log -e ./logs/02_count_reads_N4.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N4 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N18.out.log -e ./logs/02_count_reads_N18.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N18 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N19.out.log -e ./logs/02_count_reads_N19.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N19 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N20.out.log -e ./logs/02_count_reads_N20.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N20 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N21.out.log -e ./logs/02_count_reads_N21.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N21 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N25.out.log -e ./logs/02_count_reads_N25.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N25 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N26.out.log -e ./logs/02_count_reads_N26.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N26 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N27.out.log -e ./logs/02_count_reads_N27.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N27 -feature gene -attribute Name -strand_opt reverse"
