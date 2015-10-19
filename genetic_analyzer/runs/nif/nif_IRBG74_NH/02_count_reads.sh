
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Count mapped reads in gene features

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/02_count_reads_N5.out.log -e ./logs/02_count_reads_N5.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N5 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N13.out.log -e ./logs/02_count_reads_N13.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N13 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N14.out.log -e ./logs/02_count_reads_N14.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N14 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N15.out.log -e ./logs/02_count_reads_N15.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N15 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N16.out.log -e ./logs/02_count_reads_N16.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N16 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N17.out.log -e ./logs/02_count_reads_N17.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N17 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_Rhizobium_1.out.log -e ./logs/02_count_reads_Rhizobium_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples Rhizobium_1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_Rhizobium_2.out.log -e ./logs/02_count_reads_Rhizobium_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples Rhizobium_2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_Rhizobium_synnifI4_1.out.log -e ./logs/02_count_reads_Rhizobium_synnifI4_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples Rhizobium_synnifI4_1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_Rhizobium_synnifI4_2.out.log -e ./logs/02_count_reads_Rhizobium_synnifI4_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples Rhizobium_synnifI4_2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N22.out.log -e ./logs/02_count_reads_N22.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N22 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_N23.out.log -e ./logs/02_count_reads_N23.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples N23 -feature gene -attribute Name -strand_opt reverse"
