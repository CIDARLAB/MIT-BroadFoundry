
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

# A number of DE studies can be performed by varying the groups used (these are the sample column indexes from the matrix)
bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 4,5,6 -group2 2,3 -output_prefix high_vs_low_controller -bin_path $BIN_PATH/"
