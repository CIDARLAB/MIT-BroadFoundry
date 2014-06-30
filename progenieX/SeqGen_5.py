# Functions to Import from Modules
from random import random
from re_eraser import re_erase

def maine():
    seqgen(35, 40, 15, 10, 0.01, 100, 100)

def seqgen(pctA, pctT, pctC, pctG, var, seqlen, pronumber):
    
    # Define iterating list ATCG
    ATCG = ["A", "T", "C", "G"]
    
    # Define input dictionary
    pD = {"A": pctA, "T": pctT, "C": pctC, "G": pctG}

    # Define tolerance in fraction instead of percent
    tol = var / 100.
    
    # Define empty dictionaries
    dD = {}     # Binning dictionary for determining which nucleotide is substituted 
                # for random number generated below
    fD = {}     # Nucleotide fraction dictionary
    nD = {}     # Nucleotide number dictionary
    
    llimD = {}  # Lower allowable limit of variation from input nucleotide percentage
    ulimD = {}  # Upper allowable limit of variation from input nucleotide percentage

    gcD = {}    # Generated sequence nucleotide count dictionary
    gfD = {}    # Generated sequence nucleotide fraction dictionary
    vote = {}   # Check to see if sequence is within variance limits
        
    # Set promoter number to zero
    pCount = 0

    # Counter for distribution bins set to zero
    dVal = 0.

    # Populate dictionaries
    for x in ATCG:
        fD[x] = pD[x]/100. 
        nD[x] = fD[x]*seqlen
        
        llimD[x] = fD[x] - tol
        ulimD[x] = fD[x] + tol
        
        dD[x] = fD[x] + dVal
        dVal = dD[x]

    check1 = percent_error(dVal)    
    if check1 != 1 :
            print check1
            return
        
    # Open and clear output files
    pro_output = open('seqsgen.txt', 'w')
    pro_output.write('')
    pro_output.close()
    pro_output = open('seqsgen.txt', 'a')
    
    # Loop that creates number of desired sequences
    while pCount < pronumber :
        
        # Define temporary list to append nucleotides on to
        templist = []

        # Counter for variance analysis
        seq_ok = 0
        
        # Loop that generates sequences of desired length
        for nucleotide in xrange(seqlen) :

            number = random()

            if number <= dD['A'] :
                templist.append('A')
            elif dD['A'] < number <= dD['T'] :         
                templist.append('T')
            elif dD['T'] < number <= dD['C'] :         
                templist.append('C')            
            elif dD['C'] < number <= dD['G'] :         
                templist.append('G')
        
        for x in ATCG :

            gcD[x] = templist.count(x)
            gfD[x] = gcD[x]/float(seqlen)

            if llimD[x] <= gfD[x] <= ulimD[x] :
                vote[x] = 1
            else :
                vote[x] = 0

            seq_ok = vote[x] + seq_ok
        
        if seq_ok is 4 :
        
            seq = ''.join(templist)

            seq = re_erase(seq)

            pCount = pCount + 1
            pro_output.write(">EMY%(name)-spro\n%(promoter)-s\n" % {"name": pCount, "promoter": seq})
            
    pro_output.close()

def percent_error(dVal):
    # Test for proper nucleotide percentage definitions

    errmsg = "ERROR - Nucleotide percentages invalid"
    
    if abs(1-dVal) < 0.01 :
        errmsg = 1

    return errmsg

if __name__ == "__maine__" :
    maine()
