
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

# A number of DE studies can be performed by varying the groups used (these are the sample column indexes from the matrix)
bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 2 -group2 1 -output_prefix fixing_vs_nonfixing -bin_path $BIN_PATH/"
