'''
This function removes unwanted sites from a given sequence.  In PRO_Genie, it is used to remove any possible TypeIIS sites
in the promoter, ensuring efficient cloning.
'''

def maine() :
    re_erase(raw_input())
    
def re_erase(seq) :
    
    # Define unwanted restriction enzyme sites
    bpiI_f = "GAAGAC"
    bpiI_r = "GTCTTC"
    bsaI_f = "GGTCTC"
    bsaI_r = "GAGACC"
    sapI_f = "GCTCTTC"
    sapI_r = "GAAGAGC"
    mlyI_f = "GAGTC"
    mlyI_r = "GACTC"

    fix_bpiI_f = "GATGAC"
    fix_bpiI_r = "GTCATC"
    fix_bsaI_f = "GGTCAC"
    fix_bsaI_r = "GTGACC"
    fix_sapI_f = "GCTCATC"
    fix_sapI_r = "GATGAGC"
    fix_mlyI_f = "GTGTC"
    fix_mlyI_r = "GACAC"

    relist = ["bpiI_f", "bpiI_r", "bsaI_f", "bsaI_r", "sapI_f", "sapI_r", "mlyI_f", "mlyI_r"]
    reD = {"bpiI_f": bpiI_f, "bpiI_r": bpiI_r, "bsaI_f": bsaI_f, "bsaI_r":bsaI_r, "sapI_f":sapI_f, "sapI_r":sapI_r, "mlyI_f":mlyI_f, "mlyI_r":mlyI_r}
    fixD = {"bpiI_f": fix_bpiI_f, "bpiI_r": fix_bpiI_r, "bsaI_f": fix_bsaI_f, "bsaI_r":fix_bsaI_r, "sapI_f":fix_sapI_f, "sapI_r":fix_sapI_r, "mlyI_f":fix_mlyI_f, "mlyI_r":fix_mlyI_r}
    
    re_count = 0
    resD = {}
    re_seq = {}
    
    for x in relist:
        if reD[x] in seq :
            seq = seq.replace(reD[x], fixD[x])
            re_count = re_count + 1
            
    return seq

if __name__ == "__maine__" :
    maine()
