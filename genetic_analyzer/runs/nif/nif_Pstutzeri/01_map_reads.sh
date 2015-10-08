
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Map raw reads

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/01_map_reads_N8.out.log -e ./logs/01_map_reads_N8.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N8"
bsub -q forest -o ./logs/01_map_reads_N9.out.log -e ./logs/01_map_reads_N9.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N9"
