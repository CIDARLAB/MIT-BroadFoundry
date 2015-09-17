
# Generate the transcription profile for each sample

BIN_PATH=/home/unix/tgorocho/part_characterizer/bin

bsub -q forest -o ./logs/03_transcription_profile_1_1.out.log -e ./logs/03_transcription_profile_1_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 1_1"
bsub -q forest -o ./logs/03_transcription_profile_2_1.out.log -e ./logs/03_transcription_profile_2_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 2_1"
bsub -q forest -o ./logs/03_transcription_profile_3_1.out.log -e ./logs/03_transcription_profile_3_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 3_1"
bsub -q forest -o ./logs/03_transcription_profile_4_1.out.log -e ./logs/03_transcription_profile_4_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 4_1"
bsub -q forest -o ./logs/03_transcription_profile_5_1.out.log -e ./logs/03_transcription_profile_5_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 5_1"
bsub -q forest -o ./logs/03_transcription_profile_6_1.out.log -e ./logs/03_transcription_profile_6_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 6_1"
bsub -q forest -o ./logs/03_transcription_profile_7_1.out.log -e ./logs/03_transcription_profile_7_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 7_1"
bsub -q forest -o ./logs/03_transcription_profile_8_1.out.log -e ./logs/03_transcription_profile_8_1.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 8_1"

bsub -q forest -o ./logs/03_transcription_profile_1_2.out.log -e ./logs/03_transcription_profile_1_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 1_2"
bsub -q forest -o ./logs/03_transcription_profile_2_2.out.log -e ./logs/03_transcription_profile_2_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 2_2"
bsub -q forest -o ./logs/03_transcription_profile_3_2.out.log -e ./logs/03_transcription_profile_3_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 3_2"
bsub -q forest -o ./logs/03_transcription_profile_4_2.out.log -e ./logs/03_transcription_profile_4_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 4_2"
bsub -q forest -o ./logs/03_transcription_profile_5_2.out.log -e ./logs/03_transcription_profile_5_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 5_2"
bsub -q forest -o ./logs/03_transcription_profile_6_2.out.log -e ./logs/03_transcription_profile_6_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 6_2"
bsub -q forest -o ./logs/03_transcription_profile_7_2.out.log -e ./logs/03_transcription_profile_7_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 7_2"
bsub -q forest -o ./logs/03_transcription_profile_8_2.out.log -e ./logs/03_transcription_profile_8_2.err.log "python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples 8_2"
