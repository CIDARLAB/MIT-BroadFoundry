from random import random

def maine() :
    tata_sub(tbp, strength, parameterD)
    polyAT_tbp_sub(tbp, strength, parameterD)
    tss_sub(tss, strength, parameterD)
    kozak_sub(utr, strength, parameterD)
    
def tata_sub(seq, strength, parameterD) :
        
    # Since consensus TATA arise less frequently than they
    # appear in the Lubliner Very high Emax set, I need sub
    # them in at about 30% of the time when strength is 'VH'.
    
    tata_sub_list = []
    tata_sub = random()

    if strength is 'VH' :
        
        # Chance of subsituting in the TATA at a defined location in the core
        tata_yesD = {'loc_1' : parameterD['core_sub_prob']['tata'][0],
                     'loc_2' : parameterD['core_sub_prob']['tata'][1],
                     'loc_3' : parameterD['core_sub_prob']['tata'][2]}

        tata_sliceD = {'loc_1' : parameterD['core_slice']['tata'][0],
                       'loc_2' : parameterD['core_slice']['tata'][1],
                       'loc_3' : parameterD['core_slice']['tata'][2]}

        # Logic for inserting TATA at one of three locations
        if tata_sub <= tata_yesD['loc_1']:
        
            startslice = tata_sliceD['loc_1']

            if tata_sub <= tata_yesD['loc_2'] :
                startslice = tata_sliceD['loc_2']
                
            if tata_yesD['loc_2'] < tata_sub <= tata_yesD['loc_3'] :
                startslice = tata_sliceD['loc_3']
            
            endslice = startslice + 7

            tata = tata_sequence_generator()

            seqL = seq[:startslice-1]
            seqR = seq[endslice:]

            seq = seqL+tata+seqR

            tata_sub_list = [tata, startslice]
            
    return [seq, tata_sub_list]

def tata_sequence_generator() :
    
    # Consensus TATA is TATAWAWR, so need to write random generator to pick a
    # form of the consensus TATA randomly
    w = ['A','T']
    r = ['A','G']
    
    in1,in2,in3 = 0,0,0
        
    number1,number2,number3 = random(),random(),random()

    if number1 <= 0.5 :
        in1 = 1
    if number2 <= 0.5 :
        in2 = 1
    if number3 <= 0.5 :
        in3 = 1
            
    w1,w2,r1 = w[in1],w[in2],r[in3]

    tata = 'TATA%(W1)-sA%(W2)-s%(R1)-s' % {'W1': w1, 'W2' : w2, 'R1': r1}

    return tata

def polyAT_tbp_sub(seq, strength, parameterD):
    
    count = 0
    sites = parameterD['core_site']['pAT']
    
    polyAT_list = []

    while count < sites :
        
        slice_loc = parameterD['core_slice']['tbp_pAT']
        count = count + 1

        # This random generator decides whether to insert a poly dA:dT sequence, based on the probability defined
        # in the polyatD dictionary for each strength level.
        pAT_sub = random()
        
        # How likely poly dA:dT addition occurs
        pAT_yesD = parameterD['core_sub_prob']['tbp_pAT']
        
        if pAT_sub <= pAT_yesD[strength] :

            # I want to insert the poly dA:dT 15 bp into the tbp of the core.  This should be after any TATA that has been
            # inserted into the sequences.
            slice_loc_list = [slice_loc, slice_loc+13]

            pAT = tbp_pAT_seq_choose(strength, parameterD)
            
            L = seq[:slice_loc_list[0]]
            R = seq[slice_loc_list[1]:]
            
            seq = L+pAT[1]+R

            polyAT_list = [pAT[1], slice_loc_list[0]]

    return [seq, polyAT_list]

def tbp_pAT_seq_choose(strength, parameterD) :
    
    # This random generator decides which of the polydA:dT sequences to insert.  It will choose either
    # (0): polyT, (1): polyA, (2): random mix.  As polyA will likely decrease strength of the promoter
    # Or at least increase strength less than T or a mix), make polyA less likely as the strength increases.
    seq_choose = random()
    seq_choice_list = parameterD['tbp_pAT_choice_prob'][strength]
            
    # This code originally defines the poly dA:dT sequence as all T, then based on the
    # random number generated (seq_choose) decides to change the index in the list and
    # therefor the poly dA:dT sequence.
    index = 0
    if seq_choose <= seq_choice_list[0]:
        index = 2
    if seq_choice_list[0] < seq_choose <= seq_choice_list[1] :
        index = 1

    # Define the potential sequences to substitute. Poly dA:dT - 1 T, 1 A, 1 random mix
    pATD = parameterD['tbp_pAT_seq']

    return pATD[index]

def tss_sub(seq, strength, parameterD) :
    
    tss_list = []
    
    # Change how likely it is that TSS will be applied depending on strength of core
    number = random()
    prob_sub = parameterD['core_sub_prob']['tss'][strength]
    
    # Now write a chance of subsituting in the TSS at the expected TSS
    #(right at the core-utr junction)
    if number <= prob_sub :
        
        tssu = tss_upstream_chooser(strength, parameterD)
        tss_el = tss_el_chooser(strength, parameterD)

        tss = tssu+tss_el

        sliceloc = parameterD['core_slice']['tss']
    
        seq = seq[:sliceloc]+tss

        tss_list = [tss, sliceloc]

    return [seq, tss_list]
            
def tss_upstream_chooser(strength, parameterD) :
    
    tss_up_choose = random()
    choice_list = parameterD['tss_upstream_choice'][strength]
    
    index = 0
    if tss_up_choose <= choice_list[0] :
        index = 1
    if choice_list[0] < tss_up_choose <= choice_list[1] :
        index = 2
    if choice_list[1] < tss_up_choose <= choice_list[2] :
        index = 3
    
    tss_upstream = parameterD['tss_upstream'][index]
             
    return tss_upstream

def tss_el_chooser(strength, parameterD) :
             
    tss_el_choose = random()
    choice_list = parameterD['tss_el_choice'][strength]
    
    index = 0
    if tss_el_choose <= choice_list[0] :
        index = 1
    if choice_list[0] < tss_el_choose <= choice_list[1] :
        index = 2
    if choice_list[1] < tss_el_choose <= choice_list[2] :
        index = 3

    tss_el = parameterD['tss_el'][index]

    return tss_el
            
def kozak_sub(utr, strength, parameterD) :
    
    kozak_list = []
    
    number = random()
    prob_sub = parameterD['core_sub_prob']['kozak'][strength]
    
    # Now write a chance of subsituting in the Kozak at the end of the core
    if number <= prob_sub :
        
        koz_choose = random()
        choice_list = parameterD['kozak_choice'][strength]
        index = 0
        if koz_choose <= choice_list[0] :
            index = 1
        if choice_list[0] < koz_choose <= choice_list[1] :
            index = 2
        if choice_list[1] < koz_choose <= choice_list[2] :
            index = 3
            
        kozak = parameterD['kozak'][index]

        sliceloc = parameterD['core_slice']['kozak'][0]
    
        utr = utr[:sliceloc]+kozak

        if strength is 'VH' :
           double_sub = random()
           if double_sub <= parameterD['core_sub_prob']['kozak']['db']:
               sliceloc_2 = parameterD['core_slice']['kozak'][1]
               
               utr = utr[:sliceloc_2]+parameterD['kozak'][index]+kozak
              
        kozak_list = [kozak, sliceloc] 
            
    return [utr, kozak_list]

if __name__ == "__maine__" :
    maine()
