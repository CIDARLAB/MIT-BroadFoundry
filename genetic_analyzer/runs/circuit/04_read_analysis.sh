
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Analyse reads

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/04_read_analysis.out.log -e ./logs/04_read_analysis.err.log "python $BIN_PATH/read_analysis.py -settings ./data/settings.txt -bin_path $BIN_PATH/"
