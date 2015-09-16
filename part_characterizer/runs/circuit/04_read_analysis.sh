
# Analyse reads

BIN_PATH=/home/unix/tgorocho/part_characterizer/bin

bsub -q hour -W 4:0 -o ./logs/04_read_analysis.out.log -e ./logs/04_read_analysis.err.log "python $BIN_PATH/read_analysis.py -settings ./data/settings.txt"
