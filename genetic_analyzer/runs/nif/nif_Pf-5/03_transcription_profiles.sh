
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Generate the transcription profile for each sample

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/03_transcription_profile_N3.out.log -e ./logs/03_transcription_profile_N3.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N3"
bsub -q forest -o ./logs/03_transcription_profile_N4.out.log -e ./logs/03_transcription_profile_N4.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N4"
bsub -q forest -o ./logs/03_transcription_profile_N18.out.log -e ./logs/03_transcription_profile_N18.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N18"
bsub -q forest -o ./logs/03_transcription_profile_N19.out.log -e ./logs/03_transcription_profile_N19.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N19"
bsub -q forest -o ./logs/03_transcription_profile_N20.out.log -e ./logs/03_transcription_profile_N20.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N20"
bsub -q forest -o ./logs/03_transcription_profile_N21.out.log -e ./logs/03_transcription_profile_N21.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N21"
