
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Count mapped reads in gene features

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/02_count_reads_EcoliMG1655_LBWB_1.out.log -e ./logs/02_count_reads_EcoliMG1655_LBWB_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples EcoliMG1655_LBWB_1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_EcoliMG1655_LBWB_2.out.log -e ./logs/02_count_reads_EcoliMG1655_LBWB_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples EcoliMG1655_LBWB_2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_EcoliMG1655_synnifI4_1.out.log -e ./logs/02_count_reads_EcoliMG1655_synnifI4_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples EcoliMG1655_synnifI4_1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_EcoliMG1655_synnifI4_2.out.log -e ./logs/02_count_reads_EcoliMG1655_synnifI4_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples EcoliMG1655_synnifI4_2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N1.out.log -e ./logs/02_count_reads_N1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N2.out.log -e ./logs/02_count_reads_N2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N29.out.log -e ./logs/02_count_reads_N29.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N29 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N30.out.log -e ./logs/02_count_reads_N30.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N30 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N31.out.log -e ./logs/02_count_reads_N31.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N31 -feature gene -attribute Name -strand_opt reverse"
