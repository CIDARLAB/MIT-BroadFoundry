
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Generate the transcription profile for each sample

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/03_transcription_profile_N6.out.log -e ./logs/03_transcription_profile_N6.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N6"
bsub -q forest -o ./logs/03_transcription_profile_N7.out.log -e ./logs/03_transcription_profile_N7.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N7"
