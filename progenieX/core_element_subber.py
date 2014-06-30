from random import random
from xlrd import *
from excel_functions import cell
from tss_chooser import *

def maine() :
    tatasub(seq, iB, cS)
    kozaksub(core, strength, iB, cS)
    tsssub(tss, strength, iB, cS)
    polyAT_tbp_sub(tbp, strength, iB, uS)
    
def tatasub(seq, iB, cS) :

    number = random()

    # Now write a chance of subsituting in the TATA at a defined location in the core
    if number <= cell(iB, cS, 'B26') :
        
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

        startslice = int(cell(iB, cS, 'C18'))

        if number <= cell(iB, cS, 'D28') :
            startslice = int(cell(iB, cS, 'C19'))
        if cell(iB, cS, 'D28') < number <= cell(iB, cS, 'D29') :
            startslice = int(cell(iB, cS, 'C20'))
            
        endslice = startslice + 7
        
        tata = 'TATA%(W1)-sA%(W2)-s%(R1)-s' % {'W1': w1, 'W2' : w2, 'R1': r1}
        
        seqL = seq[:startslice-1]
        seqR = seq[endslice:]

        seq = seqL+tata+seqR
        
        with open('core_sub_record.txt', 'a') as rec:
            rec.write('     %(t)s\n' % {'t': tata})
            
    return seq

def kozaksub(core, strength, iB, cS) :

    x = strength

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
    
        core = core[:int(cell(iB, cS, 'C22'))]+kozak

        if x is 'VH' :
           double_sub = random()
           if double_sub <= cell(iB, cS, 'C35') :
              core = core[:int(cell(iB, cS, 'C23'))]+kozak_list[index]+kozak
              
        with open('core_sub_record.txt', 'a') as rec:
            rec.write('     %(k)s\n' % {'k': kozak})
            
    return core

def tsssub(tss, strength, iB, cS) :

    x = strength

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
    
        tss = tss[:int(cell(iB, cS, 'C21'))]+tssu+tss_el

        with open('core_sub_record.txt', 'a') as rec:
            rec.write('     %(t)s\n' % {'t':tssu+tss_el})

    return tss

def polyAT_tbp_sub(seq, strength, iB, uS):
    
    count = 0
    sites = 1

    while count < sites :
        
        slice_loc = 15
        count = count + 1

        # This random generator decides whether to insert a poly dA:dT sequence, based on the probability defined
        # in the polyatD dictionary for each strength level.
        polynum = random()
        
        # How likely poly dA:dT addition occurs
        polyatD = {1 : {'VH': float(cell(iB,uS,'B34')), 'H' : float(cell(iB,uS,'C34')),
                        'M': float(cell(iB,uS,'D34')), 'L': float(cell(iB,uS,'E34'))}}
        
        if polynum <= polyatD[1][strength] :

            # I want to insert the poly dA:dT 15 bp into the tbp of the core.  This should be after any TATA that has been
            # inserted into the sequences.
            psliceL = [slice_loc, slice_loc+13]
            
            # This random generator decides which of the polydA:dT sequences to insert.  It will choose either
            # (0): polyT, (1): polyA, (2): random mix.  As polyA will likely decrease strength of the promoter
            # Or at least increase strength less than T or a mix), make polyA less likely as the strength increases.
            seq_choose = random()
            seq_choiceD = {'VH' : [float(cell(iB,uS,'J45')), float(cell(iB,uS,'J46'))],
                           'H' : [float(cell(iB,uS,'L45')), float(cell(iB,uS,'L46'))],
                           'M' : [float(cell(iB,uS,'N45')), float(cell(iB,uS,'N46'))],
                           'L' : [float(cell(iB,uS,'P45')), float(cell(iB,uS,'P46'))]}
            
            # Define the potential sequences to substitute. Poly dA:dT - 1 T, 1 A, 1 random mix
            polyAT = [cell(iB,uS,'S47'), cell(iB,uS,'S45'), cell(iB,uS,'S46')]

            index = 0
            if seq_choose <= seq_choiceD[strength][0] :
                index = 1
            if seq_choiceD[strength][0] < seq_choose <= seq_choiceD[strength][1] :
                index = 2

            # This is the actual code for replacement. I didn't use .replace() because that searches the whole
            # string for instances of the 12 bp sequence.  There is a low chance it repeats in the UAS, but it
            # is nonzero.  Therefore, slicing seq on either side of the substitution location and then
            # concatenating ensures that only one substitution will be made at the desired location.
            L = seq[:psliceL[0]]
            R = seq[psliceL[1]:]
            
            seq = L+polyAT[index]+R

            with open('core_sub_record.txt', 'a') as rec:
                rec.write('     %(a)s\n' % {'a':polyAT[index]})

    return seq

if __name__ == "__maine__" :
    maine()
