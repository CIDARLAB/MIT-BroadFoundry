
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Count mapped reads in gene features

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/02_count_reads_N6.out.log -e ./logs/02_count_reads_N6.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N6 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N7.out.log -e ./logs/02_count_reads_N7.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N7 -feature gene -attribute Name -strand_opt reverse"
