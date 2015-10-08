
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Map raw reads

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/01_map_reads_EcoliMG1655_LBWB_1.out.log -e ./logs/01_map_reads_EcoliMG1655_LBWB_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples EcoliMG1655_LBWB_1"
bsub -q forest -o ./logs/01_map_reads_EcoliMG1655_LBWB_2.out.log -e ./logs/01_map_reads_EcoliMG1655_LBWB_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples EcoliMG1655_LBWB_2"
bsub -q forest -o ./logs/01_map_reads_EcoliMG1655_synnifI4_1.out.log -e ./logs/01_map_reads_EcoliMG1655_synnifI4_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples EcoliMG1655_synnifI4_1"
bsub -q forest -o ./logs/01_map_reads_EcoliMG1655_synnifI4_2.out.log -e ./logs/01_map_reads_EcoliMG1655_synnifI4_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples EcoliMG1655_synnifI4_2"
bsub -q forest -o ./logs/01_map_reads_N1.out.log -e ./logs/01_map_reads_N1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N1"
bsub -q forest -o ./logs/01_map_reads_N2.out.log -e ./logs/01_map_reads_N2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples N2"
