
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Map raw reads

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/01_map_reads_N5.out.log -e ./logs/01_map_reads_N5.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N5"
bsub -q forest -o ./logs/01_map_reads_N13.out.log -e ./logs/01_map_reads_N13.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N13"
bsub -q forest -o ./logs/01_map_reads_N14.out.log -e ./logs/01_map_reads_N14.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N14"
bsub -q forest -o ./logs/01_map_reads_N15.out.log -e ./logs/01_map_reads_N15.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N15"
bsub -q forest -o ./logs/01_map_reads_N16.out.log -e ./logs/01_map_reads_N16.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N16"
bsub -q forest -o ./logs/01_map_reads_N17.out.log -e ./logs/01_map_reads_N17.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N17"
