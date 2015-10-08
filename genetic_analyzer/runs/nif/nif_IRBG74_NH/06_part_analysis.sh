
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

# Perform part analysis for all promoters, terminators and ribozymes in GFF for given chromosomes
bsub -q forest -o ./logs/06_part_analysis.out.log -e ./logs/06_part_analysis.err.log "python $BIN_PATH/part_profile_analysis_T7.py -settings ./data/settings.txt"
