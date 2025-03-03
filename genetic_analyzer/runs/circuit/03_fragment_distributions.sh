
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Calculate the fragment size distributions (for profile normlisation)

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/03_fragment_distributions_1_1.out.log -e ./logs/03_fragment_distributions_1_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_1"
bsub -q forest -o ./logs/03_fragment_distributions_2_1.out.log -e ./logs/03_fragment_distributions_2_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_2"
bsub -q forest -o ./logs/03_fragment_distributions_3_1.out.log -e ./logs/03_fragment_distributions_3_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_3"
bsub -q forest -o ./logs/03_fragment_distributions_4_1.out.log -e ./logs/03_fragment_distributions_4_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_4"
bsub -q forest -o ./logs/03_fragment_distributions_5_1.out.log -e ./logs/03_fragment_distributions_5_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_5"
bsub -q forest -o ./logs/03_fragment_distributions_6_1.out.log -e ./logs/03_fragment_distributions_6_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_6"
bsub -q forest -o ./logs/03_fragment_distributions_7_1.out.log -e ./logs/03_fragment_distributions_7_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_7"
bsub -q forest -o ./logs/03_fragment_distributions_8_1.out.log -e ./logs/03_fragment_distributions_8_1.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples tube_8"

bsub -q forest -o ./logs/03_fragment_distributions_1_2.out.log -e ./logs/03_fragment_distributions_1_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_1"
bsub -q forest -o ./logs/03_fragment_distributions_2_2.out.log -e ./logs/03_fragment_distributions_2_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_2"
bsub -q forest -o ./logs/03_fragment_distributions_3_2.out.log -e ./logs/03_fragment_distributions_3_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_3"
bsub -q forest -o ./logs/03_fragment_distributions_4_2.out.log -e ./logs/03_fragment_distributions_4_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_4"
bsub -q forest -o ./logs/03_fragment_distributions_5_2.out.log -e ./logs/03_fragment_distributions_5_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_5"
bsub -q forest -o ./logs/03_fragment_distributions_6_2.out.log -e ./logs/03_fragment_distributions_6_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_6"
bsub -q forest -o ./logs/03_fragment_distributions_7_2.out.log -e ./logs/03_fragment_distributions_7_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_7"
bsub -q forest -o ./logs/03_fragment_distributions_8_2.out.log -e ./logs/03_fragment_distributions_8_2.err.log "python $BIN_PATH/fragment_distributions.py -settings ./data/settings.txt -samples flask_8"
