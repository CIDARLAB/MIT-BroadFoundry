'''
This function and associated helper functions handles the sequence generation and motif substution for the two UAS
in the promoter design.  

First, sequences are generated using SeqGen_5.  Then, based on intended strength and whether the UAS is 1 or 2, 
different transcription factor binding sites (TFBS) and poly dA:dT (pdW) motifs are substituted at defined locations.
Some of these motifs are postitively correlated with strength, while a few are negatively correlated with strength.
All motifs are drawn from literature and databases such as YEASTRACT.  

Once the substitutions are made, they are wiped of potential TypeIIS sites and the correct junctions and TypeIIS sites
added to the ends of the sequence for proper cloning.  Finally, the funciton uses the analysis package to provide
motif frequency and nucleotide content information in the form of txt files.

NOTE: This version of PRO_Genie has all parameters hard coded into the script.  PRO_GenieX has the xlrd functions 
that allow use of a spreadsheet to hold all parameters.
'''


import re
from random import random
from common_functions import *
from SeqGen_5 import seqgen
from re_eraser import re_erase

def maine() :
    uasgen(int(raw_input("Number (must be a multiple of 4):")))

def uasgen(uasnum):

    strengths = ['VH', 'H', 'M', 'L']
    ATCG = ['A', 'T', 'C', 'G']
    uas_num = [1,2]

    # These are all of the scars I will be using for all generated sequences
    # For UAS1, the left scar will be odd J and the right even J for each strength
    # For UAS2, the left scar will be A and the right odd J for each strength
    A = 'GTGC'
    B = 'AATG'
    J1 = 'TTCT'
    J2 = 'AAAC'
    J3 = 'ACTA'
    J4 = 'CTTA'
    J5 = 'CCGA'
    J6 = 'GATA'
    J7 = 'CCGA'
    J8 = 'GACC'

    uas_scar = {'VH': {1 : [J1, J2], 2 : [A, J1]},
                'H': {1 : [J3, J4], 2 : [A, J3]},
                'M': {1 : [J5, J6], 2 : [A, J5]},
                'L': {1 : [J7, J8], 2 : [A, J8]}}

    # These are the Type IIs sites to add to the end of the sequence for cloning
    bbs1_F = 'GCAGAAGACTA'
    bbs1_R = 'TAGTCTTCTGC'
            
    # Since the generator will create VH, H, L and M for each uasnum, to output the
    # number of desired sequences, the input has to be divided by 4.
    uass = uasnum/4

    for y in uas_num :

        uasD = {'VH' : [], 'H' : [], 'M' : [], 'L' : []}
                
        for x in strengths :
            # Generate random sequences with nucleotide percentages common in promoters
            seqgen(28, 34, 20, 18, 0.001, 150, uass)

            # Retrieves the generated sequences from the output file and
            # makes a list of the sequences without the names added in seqgen()
            sygenlist = lister('seqsgen.txt')
            # Eliminate the generic FASTA names applied by seqgen
            uaslist = [line for line in sygenlist if '>' not in line]
            
            # Run the elmsub() function for each sequence and eliminate TypeIIS sites that
            # may have appeared. Concatenate sequences with the appropriate scar sequence.
            for num, seq in enumerate(uaslist) :
                uas_el = elmsub(seq, x, y)
                uas = re_erase(uas_scar[x][y][0]+uas_el+uas_scar[x][y][1])
                uaslist[num] = bbs1_F+uas+bbs1_R

            # Save list to the UAS dictionary defined earlier, output in FASTA format, and analyze.
            uasD[x] = uaslist
            prefix = 'uas%(num)s' % {'num': y}
            fastaD_out(uasD, prefix)
            promoter_analysis(prefix)
    
def elmsub(seq, strength, uas) :
    
    # How many TFBS & polyAT to have a chance to add
    tfsiteD = {1 : 1, 2 : 4}
    polyATsiteD = {1 : 4, 2 : 1}

    # Where to add elements (polyAT before site & TFBS after site for UAS2
    # (polyAT after site & TFBS before site for UAS1)
    
    polyATlociD = {1 : [25, 60, 95, 130], 2 : [25]}

    seq = polyAT_sub(uas, seq, strength, polyATsiteD[uas], polyATlociD[uas])
    
    tflociD = {1 : [25], 2 : [25, 60, 95, 125]}
    
    seq = tfsub(uas, seq, strength, tfsiteD[uas], tflociD[uas])
    
    return seq

def tfsub(uas, seq, strength, sites, locilist) :

    count = 0
    x = strength
    
    while count < sites :

        slice_loc = locilist[count]
        count = count +1
        
        # This random generator decides whether to insert a TFBS sequence, based on the probability defined
        # in the tfD dictionary for each strength level.            
        tfnum = random()
        # How likley TFBS addition occurs
        tfD = {1 : {'VH': 1, 'H' : 1, 'M': 1, 'L': 0.5}, 2 : {'VH': 1, 'H' : 0.75, 'M': 0.25, 'L': 0.1}}

        if tfnum <= tfD[uas][x] :

            # Generate another random number to choose which TF to substitute, likelihood of choosing a particular
            # transcription factor is determined in the tf_probs list. (0): REB1, (0-1): RAP1, (1-2): GCR1, (2-3): ABF1,
            #(3-4): MCM1, (4): RSC3.
            tf_choose = random()    
            tf_probs = [0.25, 0.45, 0.6, 0.8, 0.95]

            if tf_choose <= tf_probs[0] :
                tf = reb1_site_chooser(x)  
            if tf_probs[0] < tf_choose <= tf_probs[1] :
                tf = rap1_site_chooser(x)
            if tf_probs[1] < tf_choose <= tf_probs[2] :
                tf = gcr1_site_chooser(x)
            if tf_probs[2] < tf_choose <= tf_probs[3] :
                tf = abf1_site_chooser(x)
            if tf_probs[3] < tf_choose <= tf_probs[4] :
                tf = mcm1_site_chooser(x)
            if tf_probs[4] < tf_choose <= 1 :
                tf = rsc3_site(x)
                
            # In UAS2, the polyAT will be before the TFBS.  In UAS1, the TFBS will be before the polyAT site.
            # The purpose of this is to wrap the TFBS region in nucleosome-disfavoring sequence and also to
            # separate the TFBS region further from the core TSS, in the manner of native promoters.
            tfsliceD = {1: [slice_loc-len(tf), slice_loc], 2 : [slice_loc, slice_loc+len(tf)]}
            
            print 'tf', count, tfsliceD[uas][0], tf
            
            # This is the actual code for replacement. I didn't use .replace() because that searches the whole
            # string for instances of the 12 bp sequence.  There is a low chance it repeats in the UAS, but it
            # is nonzero.  Therefore, slicing seq on either side of the substitution location and then
            # concatenating ensures that only one substitution will be made at the desired location.
            L = seq[:tfsliceD[uas][0]]
            R = seq[tfsliceD[uas][1]:]

            seq = L+tf+R

    return seq
            
def reb1_site_chooser(strength):
    
    # BINDING SITES FROM YEASTRACT AND MOGNO ET AL.
    reb1 = ['TTACCCGT', 'TCACCCGT', 'CAGCCCTT', 'TTACCCGG']

    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # for VH strength, only the first or second REB1 site is possible to choose.  For H strength, three are possible, adding
    # in the consensus site published in YEASTRACT.  For M, the 'weak' site from Mogno is included as a choice. For L, the
    # weak site is most likely.
    reb1_site_choiceD = {'VH' : [0.5, 0.5, 0.5], 'H' : [0.25, 0.25, 0.5], 'M' : [0.25, 0.5, 0.75], 'L' : [0.05, 0.90, 0.95]}

    reb1_choose = random()
    index = 0
    if reb1_choose <= reb1_site_choiceD[strength][0] :
        index = 1
    if reb1_site_choiceD[strength][0] < reb1_choose <= reb1_site_choiceD[strength][1] :
        index = 2
    if reb1_site_choiceD[strength][1] < reb1_choose <= reb1_site_choiceD[strength][2] :
        index = 3

    reb1_choice = reb1[index]

    return reb1_choice

def rap1_site_chooser(strength):
    
    # BINDING SITES FROM YEASTRACT AND MOGNO ET AL.
    rap1 = ['ACACCCAAGCAT', 'ACACCTGGACAT', 'ACCCCTTTTTTAC']

    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # Very similar to the code for the REB1_chooser, so check that out for more discussion.
    rap1_site_choiceD = {'VH' : [0.5, 0.5], 'H' : [0.33, 0.67], 'M' : [0.25, 0.75], 'L' : [0.15, 0.85]}

    rap1_choose = random()
    index = 0
    if rap1_choose <= rap1_site_choiceD[strength][0] :
        index = 1
    if rap1_site_choiceD[strength][0] < rap1_choose <= rap1_site_choiceD[strength][1] :
        index = 2

    rap1_choice = rap1[index]

    return rap1_choice

def gcr1_site_chooser(strength):

    # BINDING SITES FROM YEASTRACT, MOGNO ET AL., AND BAI ET AL.
    gcr1 = ['CGACTTCCT', 'CGGCATCCA', 'CAGCTTCCT', 'CAACGGAAG']
    
    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # Very similar to the code for the REB1_chooser, so check that out for more discussion.
    gcr1_site_choiceD = {'VH' : [0.33, 0.33, 0.66], 'H' : [0.25, 0.5, 0.75], 'M' : [0.2, 0.6, 0.8], 'L' : [0.1, 0.8, 0.9]}

    gcr1_choose = random()
    index = 0
    if gcr1_choose <= gcr1_site_choiceD[strength][0] :
        index = 1
    if gcr1_site_choiceD[strength][0] < gcr1_choose <= gcr1_site_choiceD[strength][1] :
        index = 2
    if gcr1_site_choiceD[strength][1] < gcr1_choose <= gcr1_site_choiceD[strength][2] :
        index = 3

    gcr1_choice = gcr1[index]

    return gcr1_choice

def abf1_site_chooser(strength):

    # BINDING SITES FROM YEASTRACT, LAST TWO BASED OFF "RTCRYYYNNNACG" 
    abf1 = ['AGCCGTAAATAGTTATCTTCCAAG', 'ATCATCTATCACG', 'GTCATTTTACACG']

    abf1_choose = random()
    index = 0
    if abf1_choose <= 0.33 :
        index = 1
    if 0.33 < abf1_choose <= 0.66 :
        index = 2
        
    abf1_choice = abf1[index]
    
    return abf1_choice

def mcm1_site_chooser(strength):

    # BINDING SITES FROM BAI ET AL.
    mcm1 = ['TTTCCGAAAACGGAAAT', 'ATACCAAATACGGTAAT']
    
    mcm1_choose = random()
    index = 0
    if mcm1_choose <= 0.5 :
        index = 1
        
    mcm1_choice = mcm1[index]
    
    return mcm1_choice

def rsc3_site(strength):

    # BINDING SITE FROM BAI ET AL.
    rsc3 = 'CGCGC'
    
    return rsc3

def polyAT_sub(uas, seq, strength, sites, locilist):
    
    count = 0

    while count < sites :
        
        slice_loc = locilist[count]
        count = count + 1

        # This random generator decides whether to insert a poly dA:dT sequence, based on the probability defined
        # in the polyatD dictionary for each strength level.
        polynum = random()
        
        # How likely poly dA:dT addition occurs
        polyatD = {1 : {'VH': 1, 'H' : 0.75, 'M': 0.25, 'L': 0.1}, 2 : {'VH': 1, 'H' : 1, 'M': 1, 'L': 1}}
        
        if polynum <= polyatD[uas][strength] :

            # I want to insert the poly dA:dT before the probable TFBS location, so the slice location chosen is
            # the endslice for the polyAT.  It will be the startslice for the TFBS.
            psliceD = {1 : [slice_loc, slice_loc+13], 2 :  [slice_loc-13, slice_loc]}
            
            # This random generator decides which of the polydA:dT sequences to insert.  It will choose either
            # (0): polyT, (1): polyA, (2): random mix A, (3): random mix B.  As polyA will likely decrease
            # strength of the promoter (or at least increase strength less than T or a mix), make polyA less likely
            # as the strength increases.
            seq_choose = random()
            seq_choiceD = {'VH' : [0., 0.2, 0.4], 'H' : [0.1, 0.25, 0.4], 'M' : [0.25, 0.5, 0.75], 'L' : [0.5, 0.75, 1]}
            
            # Define the potential sequences to substitute. Poly dA:dT - 1 T, 1 A, 2 random mix
            polyAT = ['TTTTTTTTTTTTT', 'AAAAAAAAAAAAA', 'TTAATTTAATTTT', 'ATATATTTTTAAT']

            index = 0
            if seq_choose <= seq_choiceD[strength][0] :
                index = 1
            if seq_choiceD[strength][0] < seq_choose <= seq_choiceD[strength][1] :
                index = 2
            if seq_choiceD[strength][1] < seq_choose <= seq_choiceD[strength][2] :
                index = 3

            print 'pdW', count, psliceD[uas][0], polyAT[index]
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
