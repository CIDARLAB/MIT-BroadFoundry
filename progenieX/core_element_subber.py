from random import random
from xlrd import *
from excel_functions import cell

def maine() :
    tata_sub(tbp, strength)
    polyAT_tbp_sub(tbp, strength)
    tss_sub(tss, strength)
    kozak_sub(utr, strength)
    
def tata_sub(seq, strength) :
    
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    # Since consensus TATA arise less frequently than they
    # appear in the Lubliner Very high Emax set, I need sub
    # them in at about 30% of the time when strength is 'VH'.
    
    tata_sub_list = []
    number = random()

    if strength is 'VH' :
        
        # Chance of subsituting in the TATA at a defined location in the core
        if number <= cell(iB, cS, 'B26') :
        
            startslice = int(cell(iB, cS, 'C18'))

            if number <= cell(iB, cS, 'D28') :
                startslice = int(cell(iB, cS, 'C19'))
            if cell(iB, cS, 'D28') < number <= cell(iB, cS, 'D29') :
                startslice = int(cell(iB, cS, 'C20'))
            
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

def polyAT_tbp_sub(seq, strength):
    
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    count = 0
    sites = int(cell(iB,cS,'B43'))
    
    polyAT_list = []

    while count < sites :
        
        slice_loc = int(cell(iB,cS,'C20'))
        count = count + 1

        # This random generator decides whether to insert a poly dA:dT sequence, based on the probability defined
        # in the polyatD dictionary for each strength level.
        polynum = random()
        
        # How likely poly dA:dT addition occurs
        polyatD = {1 : {'VH': float(cell(iB,cS,'B46')), 'H' : float(cell(iB,cS,'C46')),
                        'M': float(cell(iB,cS,'D46')), 'L': float(cell(iB,cS,'E46'))}}
        
        if polynum <= polyatD[1][strength] :

            # I want to insert the poly dA:dT 15 bp into the tbp of the core.  This should be after any TATA that has been
            # inserted into the sequences.
            psliceL = [slice_loc, slice_loc+13]

            # This random generator decides which of the polydA:dT sequences to insert.  It will choose either
            # (0): polyT, (1): polyA, (2): random mix.  As polyA will likely decrease strength of the promoter
            # Or at least increase strength less than T or a mix), make polyA less likely as the strength increases.
            seq_choose = random()
            seq_choiceD = {'VH' : [float(cell(iB,cS,'H40')), float(cell(iB,cS,'H41'))],
                           'H' : [float(cell(iB,cS,'J40')), float(cell(iB,cS,'J41'))],
                           'M' : [float(cell(iB,cS,'L40')), float(cell(iB,cS,'L41'))],
                           'L' : [float(cell(iB,cS,'N40')), float(cell(iB,cS,'N41'))]}
            
            # Define the potential sequences to substitute. Poly dA:dT - 1 T, 1 A, 1 random mix
            polyAT = [cell(iB,cS,'R42'), cell(iB,cS,'R41'), cell(iB,cS,'R40')]

            # This code originally defines the poly dA:dT sequence as all T, then based on the
            # random number generated (seq_choose) decides to change the index in the list and
            # therefor the poly dA:dT sequence.
            index = 0
            if seq_choose <= seq_choiceD[strength][0]:
                index = 2
            if seq_choiceD[strength][0] < seq_choose <= seq_choiceD[strength][1] :
                index = 1

            # This is the actual code for replacement. I didn't use .replace() because that searches the whole
            # string for instances of the 12 bp sequence.  There is a low chance it repeats in the UAS, but it
            # is nonzero.  Therefore, slicing seq on either side of the substitution location and then
            # concatenating ensures that only one substitution will be made at the desired location.
            L = seq[:psliceL[0]]
            R = seq[psliceL[1]:]
            
            seq = L+polyAT[index]+R

            polyAT_list = [polyAT[index], psliceL[0]]

    return [seq, polyAT_list]

def tss_sub(seq, strength) :
    
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    x = strength
    
    tss_list = []
    
    # Change how likely it is that TSS will be applied depending on strength of core
    number = random()
    prob_sub =  {'VH' : cell(iB, cS, 'C36'),
                 'H' : cell(iB, cS, 'C37'),
                 'M' : cell(iB, cS, 'C38'),
                 'L' : cell(iB, cS, 'C39')}
    
    # Now write a chance of subsituting in the TSS at the expected TSS
    #(right at the core-utr junction)
    if number <= prob_sub[x] :
        
        tssu = tss_upstream_chooser(x, iB, cS)
        tss_el = tss_el_chooser(x, iB, cS)

        tss = tssu+tss_el

        sliceloc = int(cell(iB,cS,'C21'))
    
        seq = seq[:sliceloc]+tss

        tss_list = [tss, sliceloc]

    return [seq, tss_list]
            
def tss_upstream_chooser(strength, iB, cS) :

    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS upstream has three positively correlated elements and 1 negative
    tss_upstream_list = [cell(iB,cS,'R28'), cell(iB,cS,'R29'), cell(iB,cS,'R30'), cell(iB,cS,'R31')]
        
    tssuD = {'VH' : [cell(iB,cS,'I29'), cell(iB,cS,'I30'), cell(iB,cS,'I31')],
             'H' : [cell(iB,cS,'K29'), cell(iB,cS,'K30'), cell(iB,cS,'K31')],
             'M' : [cell(iB,cS,'M29'), cell(iB,cS,'M30'), cell(iB,cS,'M31')],
             'L' : [cell(iB,cS,'O29'), cell(iB,cS,'O30'), cell(iB,cS,'O31')]}

    x = strength
    tss_up_choose = random()
    index = 0
    if tss_up_choose <= tssuD[x][0] :
        index = 1
    if tssuD[x][0] < tss_up_choose <= tssuD[x][1] :
        index = 2
    if tssuD[x][1] < tss_up_choose <= tssuD[x][2] :
        index = 3

    tss_upstream = tss_upstream_list[index]
             
    return tss_upstream

def tss_el_chooser(strength, iB, cS) :
             
    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS elements are all positively correlated, but (0) and (1) are more strongly correlated.
    # Therefore I decreased the likelihood of choosing those elements as strength decreased.
    tss_el_list = [cell(iB,cS,'R33'), cell(iB,cS,'R34'), cell(iB,cS,'R35'), cell(iB,cS,'R36')]
             
    tss_elD = {'VH' : [cell(iB,cS,'I34'), cell(iB,cS,'I35'), cell(iB,cS,'I36')],
             'H' : [cell(iB,cS,'K34'), cell(iB,cS,'K35'), cell(iB,cS,'K36')],
             'M' : [cell(iB,cS,'M34'), cell(iB,cS,'M35'), cell(iB,cS,'M36')],
             'L' : [cell(iB,cS,'O34'), cell(iB,cS,'O35'), cell(iB,cS,'O36')]}

    x = strength
    tss_el_choose = random()
    index = 0
    if tss_el_choose <= tss_elD[x][0] :
        index = 1
    if tss_elD[x][0] < tss_el_choose <= tss_elD[x][1] :
        index = 2
    if tss_elD[x][1] < tss_el_choose <= tss_elD[x][2] :
        index = 3

    tss_el = tss_el_list[index]

    return tss_el
            
def kozak_sub(core, strength) :
    
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    x = strength

    kozak_list = []

    # Change how likely it is that Kozak will be applied depending on strength of core
    number = random()
    prob_sub =  {'VH' : cell(iB, cS, 'C31'),
                 'H' : cell(iB, cS, 'C32'),
                 'M' : cell(iB, cS, 'C33'),
                 'L' : cell(iB, cS, 'C34')}
    
    # Now write a chance of subsituting in the Kozak at the end of the core
    if number <= prob_sub[x] :
        
        # Kozak mutants derived from consensus sequences in Dvir et al. In the design,
        # the -1 nucleotide will always be A because Bscar is AATG.
        # So the Kozaks are only 9 nucleotides long. The first two Kozaks
        # are positively correlated, the second two are negatively correlated.
        
        kozak_list = [cell(iB, cS, 'R19'),
                      cell(iB, cS, 'R20'),
                      cell(iB, cS, 'R21'),
                      cell(iB, cS, 'R22')]
        
        kozD = {'VH' : [cell(iB, cS, 'I20'), cell(iB, cS, 'I21'), cell(iB, cS, 'I22')],
                'H' : [cell(iB, cS, 'K20'), cell(iB, cS, 'K21'), cell(iB, cS, 'K22')],
                'M' : [cell(iB, cS, 'M20'), cell(iB, cS, 'M21'), cell(iB, cS, 'M22')],
                'L' : [cell(iB, cS, 'O20'), cell(iB, cS, 'O21'), cell(iB, cS, 'O22')]}
    
        koz_choose = random()
        index = 0
        if koz_choose <= kozD[x][0] :
            index = 1
        if kozD[x][0] < koz_choose <= kozD[x][1] :
            index = 2
        if kozD[x][1] < koz_choose <= kozD[x][2] :
            index = 3
            
        kozak = kozak_list[index]

        sliceloc = int(cell(iB, cS, 'C22'))
    
        core = core[:sliceloc]+kozak

        if x is 'VH' :
           double_sub = random()
           if double_sub <= cell(iB, cS, 'C35') :
               sliceloc_2 = int(cell(iB, cS, 'C23'))
               
               core = core[:sliceloc_2]+kozak_list[index]+kozak
              
        kozak_list = [kozak, sliceloc] 
            
    return [core, kozak_list]

if __name__ == "__maine__" :
    maine()
