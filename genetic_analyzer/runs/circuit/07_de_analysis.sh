
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

# A number of DE studies can be performed by varying the groups used (these are the sample column indexes from the matrix)
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 1,2,3,4,5,6,7,8 -output_prefix flask_vs_tube -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 13,14,15,16 -group2 9,10,11,12 -output_prefix ara_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 5,6,7,8 -group2 1,2,3,4 -output_prefix ara_comp_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 10,12,14,16 -group2 9,11,13,15 -output_prefix iptg_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 2,4,6,8 -group2 1,3,5,7 -output_prefix iptg_comp_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 11,12,15,16 -group2 9,10,13,14 -output_prefix atc_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 3,4,7,8 -group2 1,2,5,6 -output_prefix atc_comp_flask -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 14,15,16 -group2 6,7,8 -output_prefix broken_flask_vs_tube -bin_path $BIN_PATH/"

bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 10,11,12,14,15,16 -group2 9,13 -output_prefix AmtR_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 11,12,13,14,15,16 -group2 9,10 -output_prefix LitR_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 13,14,15,16 -group2 9,10,11,12 -output_prefix BM3R1_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13 -group2 14,15,16 -output_prefix SrpR_comp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,14,15,16 -group2 11,12,13 -output_prefix PhlF_comp_tube -bin_path $BIN_PATH/"


bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 9 -output_prefix state_1_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 10 -output_prefix state_2_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 11 -output_prefix state_3_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 12 -output_prefix state_4_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 13 -output_prefix state_5_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 14 -output_prefix state_6_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 15 -output_prefix state_7_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 16 -output_prefix state_8_tube -bin_path $BIN_PATH/"


bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10 -group2 9 -output_prefix rep_2_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10 -group2 10 -output_prefix rep_3_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10 -group2 11,12,13 -output_prefix rep_3_yfp_tube -bin_path $BIN_PATH/"
bsub -q forest -o ./logs/07_de_analysis.out.log -e ./logs/07_de_analysis.err.log "python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10 -group2 14,15,16 -output_prefix rep_4_tube -bin_path $BIN_PATH/"
