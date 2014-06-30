from random import random
from reverser import rev_comp

def maine() :
    re_eraser(raw_input())
    atg_eraser(seq, strength)
    nab_nrd_eraser(seq, strength)

def re_eraser(seq) :
    
    # Define unwanted restriction enzyme sites
    bpiI_f = "GAAGAC"
    bpiI_r = rev_comp(bpiI_f)
    bsaI_f = "GGTCTC"
    bsaI_r = rev_comp(bsaI_f)
    sapI_f = "GCTCTTC"
    sapI_r = rev_comp(sapI_f)
    mlyI_f = "GAGTC"
    mlyI_r = rev_comp(mlyI_f)

    fix_bpiI_f = "GATGAC"
    fix_bpiI_r = "GTCATC"
    fix_bsaI_f = "GGTCAC"
    fix_bsaI_r = "GTGACC"
    fix_sapI_f = "GCTCATC"
    fix_sapI_r = "GATGAGC"
    fix_mlyI_f = "GTGTC"
    fix_mlyI_r = "GACAC"

    relist = ["bpiI_f", "bpiI_r", "bsaI_f", "bsaI_r",
              "sapI_f", "sapI_r", "mlyI_f", "mlyI_r"]
    
    reD = {"bpiI_f": bpiI_f, "bpiI_r": bpiI_r,
           "bsaI_f": bsaI_f, "bsaI_r":bsaI_r,
           "sapI_f":sapI_f, "sapI_r":sapI_r,
           "mlyI_f":mlyI_f, "mlyI_r":mlyI_r}
    
    fixD = {"bpiI_f": fix_bpiI_f, "bpiI_r": fix_bpiI_r,
            "bsaI_f": fix_bsaI_f, "bsaI_r":fix_bsaI_r,
            "sapI_f":fix_sapI_f, "sapI_r":fix_sapI_r,
            "mlyI_f":fix_mlyI_f, "mlyI_r":fix_mlyI_r}
    
    re_count = 0
    resD = {}
    re_seq = {}
    
    for x in relist:
        if reD[x] in seq :
            seq = seq.replace(reD[x], fixD[x])
            re_count = re_count + 1
            
    return seq

def atg_eraser(seq, strength) :

    s = strength
    
    atg = "ATG"

    if s is not "L":

        fix = "AAG"
    
        x = random()

        if 0 <= x <= 0.25:
            fix = "TTG"
        if 0.25 <= x <= 0.5:
            fix = "ATC"
        
        if atg in seq :
            seq = seq.replace(atg, fix)
            
    return seq

def nab_nrd_eraser(seq, strength) :

    s = strength
    
    nab3 = "TCTT"
    nrd1_1 = "GTAA"
    nrd1_2 = "GTAG"

    if s is not "L":

        nab_fix = "ACTA"
        nrd_fix = "GTAT"
        
        if nab3 in seq :
            seq = seq.replace(nab3, nab_fix)
        if nrd1_1 in seq :
            seq = seq.replace(nrd1_1, nrd_fix)
        if nrd1_2 in seq :
            seq = seq.replace(nrd1_2, nrd_fix)
            
    return seq

if __name__ == "__maine__" :
    maine()
