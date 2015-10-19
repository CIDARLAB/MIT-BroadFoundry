
#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Generate the transcription profile for each sample

BIN_PATH=/home/unix/tgorocho/genetic_analyzer/bin

bsub -q forest -o ./logs/03_transcription_profile_N5.out.log -e ./logs/03_transcription_profile_N5.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N5"
bsub -q forest -o ./logs/03_transcription_profile_N13.out.log -e ./logs/03_transcription_profile_N13.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N13"
bsub -q forest -o ./logs/03_transcription_profile_N14.out.log -e ./logs/03_transcription_profile_N14.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N14"
bsub -q forest -o ./logs/03_transcription_profile_N15.out.log -e ./logs/03_transcription_profile_N15.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N15"
bsub -q forest -o ./logs/03_transcription_profile_N16.out.log -e ./logs/03_transcription_profile_N16.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N16"
bsub -q forest -o ./logs/03_transcription_profile_N17.out.log -e ./logs/03_transcription_profile_N17.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N17"
bsub -q forest -o ./logs/03_transcription_profile_Rhizobium_1.out.log -e ./logs/03_transcription_profile_Rhizobium_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples Rhizobium_1"
bsub -q forest -o ./logs/03_transcription_profile_Rhizobium_2.out.log -e ./logs/03_transcription_profile_Rhizobium_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples Rhizobium_2"
bsub -q forest -o ./logs/03_transcription_profile_Rhizobium_synnifI4_1.out.log -e ./logs/03_transcription_profile_Rhizobium_synnifI4_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples Rhizobium_synnifI4_1"
bsub -q forest -o ./logs/03_transcription_profile_Rhizobium_synnifI4_2.out.log -e ./logs/03_transcription_profile_Rhizobium_synnifI4_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples Rhizobium_synnifI4_2"
bsub -q forest -o ./logs/03_transcription_profile_N22.out.log -e ./logs/03_transcription_profile_N22.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N22"
bsub -q forest -o ./logs/03_transcription_profile_N23.out.log -e ./logs/03_transcription_profile_N23.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples N23"
