
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Map raw reads

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/01_map_reads_N3.out.log -e ./logs/01_map_reads_N3.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N3"
bsub -q forest -o ./logs/01_map_reads_N4.out.log -e ./logs/01_map_reads_N4.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N4"
bsub -q forest -o ./logs/01_map_reads_N18.out.log -e ./logs/01_map_reads_N18.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N18"
bsub -q forest -o ./logs/01_map_reads_N19.out.log -e ./logs/01_map_reads_N19.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N19"
bsub -q forest -o ./logs/01_map_reads_N20.out.log -e ./logs/01_map_reads_N20.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N20"
bsub -q forest -o ./logs/01_map_reads_N21.out.log -e ./logs/01_map_reads_N21.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N21"
