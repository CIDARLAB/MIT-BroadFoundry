
# Map raw reads

BIN_PATH=/home/unix/tgorocho/part_characterizer/bin

bsub -q forrest -o ./logs/01_map_reads_1_1.out.log -e ./logs/01_map_reads_1_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 1_1"
bsub -q forrest -o ./logs/01_map_reads_2_1.out.log -e ./logs/01_map_reads_2_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 2_1"
bsub -q forrest -o ./logs/01_map_reads_3_1.out.log -e ./logs/01_map_reads_3_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 3_1"
bsub -q forrest -o ./logs/01_map_reads_4_1.out.log -e ./logs/01_map_reads_4_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 4_1"
bsub -q forrest -o ./logs/01_map_reads_5_1.out.log -e ./logs/01_map_reads_5_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 5_1"
bsub -q forrest -o ./logs/01_map_reads_6_1.out.log -e ./logs/01_map_reads_6_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 6_1"
bsub -q forrest -o ./logs/01_map_reads_7_1.out.log -e ./logs/01_map_reads_7_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 7_1"
bsub -q forrest -o ./logs/01_map_reads_8_1.out.log -e ./logs/01_map_reads_8_1.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 8_1"

bsub -q forrest -o ./logs/01_map_reads_1_2.out.log -e ./logs/01_map_reads_1_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 1_2"
bsub -q forrest -o ./logs/01_map_reads_2_2.out.log -e ./logs/01_map_reads_2_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 2_2"
bsub -q forrest -o ./logs/01_map_reads_3_2.out.log -e ./logs/01_map_reads_3_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 3_2"
bsub -q forrest -o ./logs/01_map_reads_4_2.out.log -e ./logs/01_map_reads_4_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 4_2"
bsub -q forrest -o ./logs/01_map_reads_5_2.out.log -e ./logs/01_map_reads_5_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 5_2"
bsub -q forrest -o ./logs/01_map_reads_6_2.out.log -e ./logs/01_map_reads_6_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 6_2"
bsub -q forrest -o ./logs/01_map_reads_7_2.out.log -e ./logs/01_map_reads_7_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 7_2"
bsub -q forrest -o ./logs/01_map_reads_8_2.out.log -e ./logs/01_map_reads_8_2.err.log "python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples 8_2"
