
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

# A number of DE studies can be performed by varying the groups used (these are the sample column indexes from the matrix)
bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 1,2,3,4,5,6,7,8 -group2 9,10,11,12,13,14,15,16 -output_prefix tube_vs_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12 -group2 13,14,15,16 -output_prefix ara_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 1,2,3,4 -group2 5,6,7,8 -output_prefix ara_comp_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,11,13,15 -group2 10,12,14,16 -output_prefix iptg_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 1,3,5,7 -group2 2,4,6,8 -output_prefix iptg_comp_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,13,14, -group2 11,12,15,16 -output_prefix atc_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 1,2,5,6 -group2 3,4,7,8 -output_prefix atc_comp_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/05_de_analysis.out.log -e ./logs/05_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 6,7,8, -group2 14,15,16 -output_prefix broken_flask_vs_tube -bin_path $BIN_PATH/"
