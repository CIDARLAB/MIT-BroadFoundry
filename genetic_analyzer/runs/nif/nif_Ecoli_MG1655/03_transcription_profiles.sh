
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Generate the transcription profile for each sample

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/03_transcription_profile_EcoliMG1655_LBWB_1.out.log -e ./logs/03_transcription_profile_EcoliMG1655_LBWB_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples EcoliMG1655_LBWB_1"
bsub -q forest -o ./logs/03_transcription_profile_EcoliMG1655_LBWB_2.out.log -e ./logs/03_transcription_profile_EcoliMG1655_LBWB_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples EcoliMG1655_LBWB_2"
bsub -q forest -o ./logs/03_transcription_profile_EcoliMG1655_synnifI4_1.out.log -e ./logs/03_transcription_profile_EcoliMG1655_synnifI4_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples EcoliMG1655_synnifI4_1"
bsub -q forest -o ./logs/03_transcription_profile_EcoliMG1655_synnifI4_2.out.log -e ./logs/03_transcription_profile_EcoliMG1655_synnifI4_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples EcoliMG1655_synnifI4_2"
bsub -q forest -o ./logs/03_transcription_profile_N1.out.log -e ./logs/03_transcription_profile_N1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N1"
bsub -q forest -o ./logs/03_transcription_profile_N2.out.log -e ./logs/03_transcription_profile_N2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N2"
bsub -q forest -o ./logs/03_transcription_profile_N29.out.log -e ./logs/03_transcription_profile_N29.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N29"
bsub -q forest -o ./logs/03_transcription_profile_N30.out.log -e ./logs/03_transcription_profile_N30.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N30"
bsub -q forest -o ./logs/03_transcription_profile_N31.out.log -e ./logs/03_transcription_profile_N31.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N31"
