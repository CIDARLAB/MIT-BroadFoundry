
# Generate the transcription profile for each sample

bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 1_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 2_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 3_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 4_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 5_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 6_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 7_1"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 8_1"

bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 1_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 2_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 3_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 4_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 5_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 6_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 7_2"
bsub -q hour -W 4:0 -N "python transcription_profile.py -settings ./data/settings.txt -samples 8_2"
