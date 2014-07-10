from random import random
from xlrd import *
from common_functions import *
from excel_functions import cell
from tfbs_selector import *

def maine() :
    elmsub(seq, strength, uas, iB, uS)
    tfsub(uas, seq, strength, sites, locilist)
    polyAT_sub(uas, seq, strength, sites, locilist)
    
def elmsub(seq, strength, uas, iB, uS) :
    
    # How many TFBS & polyAT to have a chance to add
    tfsiteD = {1 : int(cell(iB,uS,'A8')),
               2 : int(cell(iB,uS,'C8'))}
    
    polyATsiteD = {1 : int(cell(iB,uS,'B8')),
                   2 : int(cell(iB,uS,'D8'))}

    # Where to add elements (polyAT before site & TFBS after site for UAS2
    # (polyAT after site & TFBS before site for UAS1)

    tflociD = {1 : [int(cell(iB,uS,'A13')), int(cell(iB,uS,'A14'))],
               2 : [int(cell(iB,uS,'C13')), int(cell(iB,uS,'C14')),
                    int(cell(iB,uS,'C15')), int(cell(iB,uS,'C16'))]}
    
    seq = tfsub(uas, seq, strength, tfsiteD[uas], tflociD[uas], iB, uS)
    
    polyATlociD = {1 : [int(cell(iB,uS,'B13')), int(cell(iB,uS,'B14')),
                        int(cell(iB,uS,'B15')), int(cell(iB,uS,'B16'))],
                   2 : [int(cell(iB,uS,'D13'))]}

    seq = polyAT_sub(uas, seq, strength, polyATsiteD[uas], polyATlociD[uas], iB, uS)
    
    return seq

def tfsub(uas, seq, strength, sites, locilist, iB, uS) :

    count = 0
    x = strength
    
    while count < sites :

        slice_loc = locilist[count]
        count = count +1
        
        # This random generator decides whether to insert a TFBS sequence, based on
        #the probability defined in the tfD dictionary for each strength level.            
        tfnum = random()
        # How likley TFBS addition occurs
        tfD = {1 : {'VH': float(cell(iB,uS,'B20')), 'H' : float(cell(iB,uS,'C20')),
                    'M': float(cell(iB,uS,'D20')), 'L': float(cell(iB,uS,'E20'))},
               2 : {'VH': float(cell(iB,uS,'B21')), 'H' : float(cell(iB,uS,'C21')),
                    'M': float(cell(iB,uS,'D21')), 'L': float(cell(iB,uS,'E21'))}}
        
        # This ensures the first site is always inserted with a TF site.  For UAS2, the other
        # three are inserted with the probability in the ProGenie_Parameters.xlsx.  For UAS1,
        # there is only one other site that has the probability of insertion.
        if count =< 1 :

            # Generate another random number to choose which TF to substitute, likelihood
            #of choosing a particular transcription factor is determined in the tf_probs list.
            #(0): REB1, (0-1): RAP1, (1-2): GCR1, (2-3): ABF1,(3-4): MCM1, (4): RSC3.
            # In UAS2, to limit the number of patterns, as well as eliminate those TFs that act
            # close to the transcription start site, only allow REB1, RAP1, and GCR1 sites.
            tf_choose = random()
            
            tf_probs = {1: [float(cell(iB,uS,'C25')), float(cell(iB,uS,'C26')),
                        float(cell(iB,uS,'C27')), float(cell(iB,uS,'C28')),
                        float(cell(iB,uS,'C29'))],
                        2: [float(cell(iB,uS,'E25')), float(cell(iB,uS,'E26')),
                        float(cell(iB,uS,'E27')), 1, 1]}

            if tf_choose <= tf_probs[uas][0] :
                tf = reb1_site_chooser(x, iB, uS)
                tfn = 'REB1'
            if tf_probs[uas][0] < tf_choose <= tf_probs[uas][1] :
                tf = rap1_site_chooser(x, iB, uS)
                tfn = 'RAP1'
            if tf_probs[uas][1] < tf_choose <= tf_probs[uas][2] :
                tf = gcr1_site_chooser(x, iB, uS)
                tfn = 'GCR1'
            if tf_probs[uas][2] < tf_choose <= tf_probs[uas][3] :
                tf = abf1_site_chooser(x, iB, uS)
                tfn = 'ABF1'
            if tf_probs[uas][3] < tf_choose <= tf_probs[uas][4] :
                tf = mcm1_site_chooser(x, iB, uS)
                tfn = 'MCM1'
            if tf_probs[uas][4] < tf_choose <= 1 :
                tf = rsc3_site(x, iB, uS)
                tfn = 'RSC3'
                
            # In UAS2, the polyAT will be before the TFBS.  In UAS1, the TFBS will be before the polyAT site.
            # The purpose of this is to wrap the TFBS region in nucleosome-disfavoring sequence and also to
            # separate the TFBS region further from the core TSS, in the manner of native promoters.
            tfsliceD = {1: [slice_loc-len(tf), slice_loc], 2 : [slice_loc, slice_loc+len(tf)]}

            with open('uas%(y)s_sub_record.txt' % {'y': uas}, 'a') as rec :
                rec.write('tf %(c)s %(l)s %(t)s %(s)s\n' % {'c':count, 'l':tfsliceD[uas][0], 't': tfn, 's':tf})
            
            # This is the actual code for replacement. I didn't use .replace() because that searches the whole
            # string for instances of the 12 bp sequence.  There is a low chance it repeats in the UAS, but it
            # is nonzero.  Therefore, slicing seq on either side of the substitution location and then
            # concatenating ensures that only one substitution will be made at the desired location.
            L = seq[:tfsliceD[uas][0]]
            R = seq[tfsliceD[uas][1]:]

            seq = L+tf+R
            
        if count > 1:
            if tfnum <= tfD[uas][x] :

                # Generate another random number to choose which TF to substitute, likelihood
                #of choosing a particular transcription factor is determined in the tf_probs list.
                #(0): REB1, (0-1): RAP1, (1-2): GCR1, (2-3): ABF1,(3-4): MCM1, (4): RSC3.
                # In UAS2, to limit the number of patterns, as well as eliminate those TFs that act
                # close to the transcription start site, only allow REB1, RAP1, and GCR1 sites.
                tf_choose = random()
            
                tf_probs = {1: [float(cell(iB,uS,'C25')), float(cell(iB,uS,'C26')),
                            float(cell(iB,uS,'C27')), float(cell(iB,uS,'C28')),
                            float(cell(iB,uS,'C29'))],
                            2: [float(cell(iB,uS,'E25')), float(cell(iB,uS,'E26')),
                            float(cell(iB,uS,'E27')), 1, 1]}

                if tf_choose <= tf_probs[uas][0] :
                    tf = reb1_site_chooser(x, iB, uS)
                    tfn = 'REB1'
                if tf_probs[uas][0] < tf_choose <= tf_probs[uas][1] :
                    tf = rap1_site_chooser(x, iB, uS)
                    tfn = 'RAP1'
                if tf_probs[uas][1] < tf_choose <= tf_probs[uas][2] :
                    tf = gcr1_site_chooser(x, iB, uS)
                    tfn = 'GCR1'
                if tf_probs[uas][2] < tf_choose <= tf_probs[uas][3] :
                    tf = abf1_site_chooser(x, iB, uS)
                    tfn = 'ABF1'
                if tf_probs[uas][3] < tf_choose <= tf_probs[uas][4] :
                    tf = mcm1_site_chooser(x, iB, uS)
                    tfn = 'MCM1'
                if tf_probs[uas][4] < tf_choose <= 1 :
                    tf = rsc3_site(x, iB, uS)
                    tfn = 'RSC3'
                
                # In UAS2, the polyAT will be before the TFBS.  In UAS1, the TFBS will be before the polyAT site.
                # The purpose of this is to wrap the TFBS region in nucleosome-disfavoring sequence and also to
                # separate the TFBS region further from the core TSS, in the manner of native promoters.
                    tfsliceD = {1: [slice_loc-len(tf), slice_loc], 2 : [slice_loc, slice_loc+len(tf)]}

                with open('uas%(y)s_sub_record.txt' % {'y': uas}, 'a') as rec :
                    rec.write('tf %(c)s %(l)s %(t)s %(s)s\n' % {'c':count, 'l':tfsliceD[uas][0], 't': tfn, 's':tf})
                
                # This is the actual code for replacement. I didn't use .replace() because that searches the whole
                # string for instances of the 12 bp sequence.  There is a low chance it repeats in the UAS, but it
                # is nonzero.  Therefore, slicing seq on either side of the substitution location and then
                # concatenating ensures that only one substitution will be made at the desired location.
                L = seq[:tfsliceD[uas][0]]
                R = seq[tfsliceD[uas][1]:]

                seq = L+tf+R

    return seq

def polyAT_sub(uas, seq, strength, sites, locilist, iB, uS):
    
    count = 0

    while count < sites :
        
        slice_loc = locilist[count]
        count = count + 1

        # This random generator decides whether to insert a poly dA:dT sequence, based on the probability defined
        # in the polyatD dictionary for each strength level.
        polynum = random()
        
        # How likely poly dA:dT addition occurs
        polyatD = {1 : {'VH': float(cell(iB,uS,'B34')), 'H' : float(cell(iB,uS,'C34')),
                        'M': float(cell(iB,uS,'D34')), 'L': float(cell(iB,uS,'E34'))},
                   2 : {'VH': float(cell(iB,uS,'B35')), 'H' : float(cell(iB,uS,'C35')),
                        'M': float(cell(iB,uS,'D35')), 'L': float(cell(iB,uS,'E35'))}}
        
        if polynum <= polyatD[uas][strength] :

            # I want to insert the poly dA:dT before the probable TFBS location, so the slice location chosen is
            # the endslice for the polyAT.  It will be the startslice for the TFBS.
            psliceD = {1 : [slice_loc, slice_loc+13], 2 :  [slice_loc-13, slice_loc]}
            
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

            with open('uas%(y)s_sub_record.txt' % {'y': uas}, 'a') as rec :
                rec.write('pdW %(c)s %(l)s %(s)s\n' % {'c':count, 'l':psliceD[uas][0], 's':polyAT[index]})
            
            # This is the actual code for replacement. I didn't use .replace() because that searches the whole
            # string for instances of the 12 bp sequence.  There is a low chance it repeats in the UAS, but it
            # is nonzero.  Therefore, slicing seq on either side of the substitution location and then
            # concatenating ensures that only one substitution will be made at the desired location.
            L = seq[:psliceD[uas][0]]
            R = seq[psliceD[uas][1]:]
            
            seq = L+polyAT[index]+R

    return seq

if __name__ == "__maine__" :
    maine()
