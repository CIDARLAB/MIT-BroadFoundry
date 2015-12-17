
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

# Perform part analysis for all promoters, terminators and ribozymes in GFF for given chromosomes
bsub -q forest -o ./logs/09_promoter_fitting_tubes.out.log -e ./logs/09_promoter_fitting_tubes.err.log "python $BIN_PATH/promoter_fitting.py -settings ./data/settings.txt -output_name tube -samples tube_1,tube_2,tube_3,tube_4,tube_5,tube_6,tube_7,tube_8"
bsub -q forest -o ./logs/09_promoter_fitting_flasks.out.log -e ./logs/09_promoter_fitting_flasks.err.log "python $BIN_PATH/promoter_fitting.py -settings ./data/settings.txt -output_name flask -samples flask_1,flask_2,flask_3,flask_4,flask_5,flask_6,flask_7,flask_8"
