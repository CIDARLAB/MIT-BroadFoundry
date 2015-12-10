
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Generate the transcription profile for each sample

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/04_transcription_profile_1_1.out.log -e ./logs/04_transcription_profile_1_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_1"
bsub -q forest -o ./logs/04_transcription_profile_2_1.out.log -e ./logs/04_transcription_profile_2_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_2"

bsub -q forest -o ./logs/04_transcription_profile_1_2.out.log -e ./logs/04_transcription_profile_1_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_1"
bsub -q forest -o ./logs/04_transcription_profile_2_2.out.log -e ./logs/04_transcription_profile_2_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_2"
