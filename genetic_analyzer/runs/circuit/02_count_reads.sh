
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Count mapped reads in gene features

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/02_count_reads_1_1.out.log -e ./logs/02_count_reads_1_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_2_1.out.log -e ./logs/02_count_reads_2_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_3_1.out.log -e ./logs/02_count_reads_3_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_3 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_4_1.out.log -e ./logs/02_count_reads_4_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_4 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_5_1.out.log -e ./logs/02_count_reads_5_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_5 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_6_1.out.log -e ./logs/02_count_reads_6_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_6 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_7_1.out.log -e ./logs/02_count_reads_7_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_7 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_8_1.out.log -e ./logs/02_count_reads_8_1.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_8 -feature gene -attribute Name -strand_opt reverse"

bsub -q forest -o ./logs/02_count_reads_1_2.out.log -e ./logs/02_count_reads_1_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_1 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_2_2.out.log -e ./logs/02_count_reads_2_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_2 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_3_2.out.log -e ./logs/02_count_reads_3_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_3 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_4_2.out.log -e ./logs/02_count_reads_4_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_4 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_5_2.out.log -e ./logs/02_count_reads_5_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_5 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_6_2.out.log -e ./logs/02_count_reads_6_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_6 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_7_2.out.log -e ./logs/02_count_reads_7_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_7 -feature gene -attribute Name -strand_opt reverse"
bsub -q forest -o ./logs/02_count_reads_8_2.out.log -e ./logs/02_count_reads_8_2.err.log "python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_8 -feature gene -attribute Name -strand_opt reverse"