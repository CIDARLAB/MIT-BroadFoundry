'''
This function searches for a concensus TATA (TATAWAWR) within sequences. It is used in the promoter analysis
group of functions.
'''

import re

def maine() :
    tatafind(nom, seq, count, step)

def tatafind(nom, seq, count, step) :

    yesF, yesR = 0, 0
    locationF, locationR = 'NA', 'NA'
    
    match_F = re.search('(TATA[AT]A[AT][AG])', seq)
    match_R = re.search('([CT][AT]T[AT]TATA)', seq)
    
    if match_F:
        yesF = 1
        loc_loc = match_F.start()
        locationF = loc_loc+count*step

    if match_R:
        yesR = 1
        loc_loc = match_R.start()
        locationR = loc_loc+count*step

    return nom, yesF, locationF, yesR, locationR

if __name__ == "__maine__" :
    maine()
