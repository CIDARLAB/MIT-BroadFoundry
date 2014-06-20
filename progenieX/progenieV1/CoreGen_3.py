'''
This is the main function for generating the 150 bp core promoters in PRO_genie.  This function is designed to 
generate three different 50 bp subsections of the core: the TATA binding portion (TBP), the probable transcription start
site (TSS), and the probable 5' UTR (UTR).  These subsections each have different nucleotide content, which also varies
with intended strength, as well as unique motifs that can be substituted.  CoreGen_3 invokes a sequence generator funciton
in SeqGen_5, then passes those sequences to motif substitution  and eraser functions.  These functions sub in TATA boxes, 
transcription start site elements, and consensus Kozak sequences, while removing undesired ATG codons from the expected 
5' UTR and removing NAB1 and NRD1 sites, proteins that have been implicated in initiating the decay of nonsense RNAS.

Once substitutions are made and the unwanted sites removed, sequences are returned and then formatted with TypeIIS sites, 
junctions for cloning, and any TypeIIS sites that may have arisen in the generation and substitution stages are removed.
The script then runs the analysis package I have written to analyze promoter sequences.

NOTE:  In this version, the parameters for the functions are hard coded into the script.  In PRO_GenieX, these values are
contained in an Excel spreadsheet that makes editing parameters much more straightforward.
'''

import re
from random import random
from common_functions import *
from SeqGen_5 import seqgen
from re_eraser import re_erase
from pro_analysis import promoter_analysis

def maine() :
    coregen(int(raw_input("Number (must be a multiple of 4):")))

def coregen(cornum) :
   
    # Define strength level dictionaries and associated fractions from Lubliner et al. 2013
    strengths = ['VH', 'H', 'M', 'L']
    subpart = ['tbp', 'tss', 'utr']
    ATCG = ['A', 'T', 'C', 'G']
    

    syxD = {'VH': {'tbp': [], 'tss': [], 'utr': []},
            'H': {'tbp': [], 'tss': [], 'utr': []},
            'M': {'tbp': [], 'tss': [], 'utr': []},
            'L': {'tbp': [], 'tss': [], 'utr': []}}

    coreD = {'VH': [], 'H': [], 'M': [], 'L': []}
    
    specD = {'VH': {'tbp': {'A': 30, 'T': 34, 'C': 18, 'G': 18},
                    'tss': {'A': 24, 'T': 48, 'C': 18, 'G': 10},
                    'utr': {'A': 40, 'T': 24, 'C': 20, 'G': 16}},
             'H': {'tbp': {'A': 32, 'T': 36, 'C': 16, 'G': 16},
                    'tss': {'A': 32, 'T': 38, 'C': 16, 'G': 14},
                    'utr': {'A': 44, 'T': 22, 'C': 18, 'G': 16}},
             'M': {'tbp': {'A': 36, 'T': 30, 'C': 16, 'G': 18},
                    'tss': {'A': 34, 'T': 30, 'C': 18, 'G': 18},
                    'utr': {'A': 36, 'T': 28, 'C': 18, 'G': 18}},
             'L': {'tbp': {'A': 34, 'T': 30, 'C': 18, 'G': 18},
                    'tss': {'A': 36, 'T': 28, 'C': 18, 'G': 18},
                    'utr': {'A': 30, 'T': 34, 'C': 18, 'G': 18}}}
    
    le = 50
    tol = 0.0001

    # These are all of the scars I will be using for all generated sequences
    # For core promoters, the left scar will be even J and right will be B
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

    core_scar = {'VH': [J2, B],
                 'H': [J4, B],
                 'M': [J6, B],
                 'L': [J8, B]}

    # These are the Type IIs sites to add to the end of the sequence for cloning
    bbs1_F = 'GCAGAAGACTA'
    bbs1_R = 'TAGTCTTCTGC'
    
    # Since the generator will create VH, H, L and M for each cornum, to output the
    # number of desired sequences, the input has to be divided by 4.
    cors = cornum/4
    
    for x in strengths : 

        for y in subpart :

            # Generate 50 bp sequences for each region
            seqgen(specD[x][y]['A'], specD[x][y]['T'], specD[x][y]['C'], specD[x][y]['G'], tol, le, cors)

            # Retrieves the generated sequences from the output file and
            # makes a list of the sequences without the names added in seqgen()
            sygenlist = lister('seqsgen.txt')
            # Eliminate the generic FASTA names applied by seqgen
            sylist = [line for line in sygenlist if '>' not in line]
            
            # Since consensus TATA arise less frequently than they appear in the Lubliner
            # Very high Emax set, I needed to write a script that subs them in at about 30%
            # of the time the tatagen() function is executed.
            # Since the consensus Kozak is also unlikely to appear at the end of the UTR,
            # I also wrote a function that substitutes it in.
            if x is 'VH' :
                if y is 'tbp':
                    for num, sy in enumerate(sylist) :
                        sy = tatasub(sy)
                        sylist[num] = sy
                        
            # Since the consensus Kozak is also unlikely to appear at the end of the UTR,
            # I also wrote a function that substitutes it in.
            if y is 'utr':
                for num, sy in enumerate(sylist) :
                    sy = kozaksub(sy, x)
                    sy = atg_eraser(sy, x)
                    sy = nab_nrd_eraser(sy, x)
                    sylist[num] = sy

            # Since there is data available about elements around the TSS from Lubliner
            # I also wrote a function that puts in TSS elements right at the end of the
            # tss element
            if y is 'tss' :
                for num, sy in enumerate(sylist) :
                    sy = tsssub(sy, x)
                    sy = nab_nrd_eraser(sy, x)
                    sylist[num] = sy

            # This populates the different segments into one dictionary
            syxD[x][y] = sylist

    # This stitches together the 50 bp segments into a full core promoter and adds TypeIIS sites.
    for x in strengths :

        coreD[x] = syxD[x]['tbp']
        
        for num, seq in enumerate(syxD[x][y]):
            core = re_erase(core_scar[x][0]+syxD[x]['tbp'][num]+syxD[x]['tss'][num]+syxD[x]['utr'][num]+core_scar[x][1])
            coreD[x][num] = bbs1_F+core+bbs1_R

    # Name output data files, format sequence dictionary into FASTA, and analyze for motifs.
    prefix = 'core'
    
    fastaD_out(coreD, prefix)    

    promoter_analysis(prefix)

def tatasub(seq) :

    number = random()

    # Now write a chance of subsituting in the TATA at a defined location in the core
    if number <= 0.25 :
        
        # Consensus TATA is TATAWAWR, so need to write random generator to pick a form of the consensus TATA randomly
        w = ['A','T']
        r = ['A','G']
    
        in1 = 0
        in2 = 0
        in3 = 0
        
        number1 = random()
        number2 = random()
        number3 = random()

        if number1 <= 0.5 :
            in1 = 1
        if number2 <= 0.5 :
            in2 = 1
        if number3 <= 0.5 :
            in3 = 1
            
        w1 = w[in1]
        w2 = w[in2]
        r1 = r[in3]

        startslice = 2

        if number <= 0.0675 :
            startslice = 7
        if 0.0675 < number <= 0.125 :
            startslice = 21
            
        endslice = startslice + 7
        
        tata = 'TATA%(W1)-sA%(W2)-s%(R1)-s' % {'W1': w1, 'W2' : w2, 'R1': r1}
        
        seqL = seq[:startslice-1]
        seqR = seq[endslice:]

        seq = seqL+tata+seqR

    return seq

def kozaksub(core, strength) :

    x = strength

    # Change how likely it is that Kozak will be applied depending on strength of core
    number = random()
    prob_sub =  {'VH' : 0.95, 'H' : 0.5, 'M' : 0.25, 'L' : 0.67}
    
    # Now write a chance of subsituting in the Kozak at the end of the core
    if number <= prob_sub[x] :
        
        # Kozak mutants derived from consensus sequences in Dvir et al. In the design, the -1 nucleotide will
        # always be A because Bscar is AATG.  So the Kozaks are only 9 nucleotides long.
        # The first two Kozaks are positively correlated, the second two are negatively correlated.
        kozak_list = ['AAAAGTAAA', 'AAAAACAAA', 'CCACCGGCG', 'CCACCAGTG'] 
        kozD = {'VH' : [0.5, 0.5, 0.5], 'H' : [0.5, 0.5, 0.5], 'M' : [0.4, 0.5, 0.6], 'L' : [0.1, 0.5, 0.9]}
    
        koz_choose = random()
        index = 0
        if koz_choose <= kozD[x][0] :
            index = 1
        if kozD[x][0] < koz_choose <= kozD[x][1] :
            index = 2
        if kozD[x][1] < koz_choose <= kozD[x][2] :
            index = 3
            
        kozak = kozak_list[index]
    
        core = core[:41]+kozak

        if x is 'VH' :
            double_sub = random()
            if double_sub <= 0.75 :
                core = core[:32]+kozak_list[1]+kozak
    
    return core

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

def tsssub(tss, strength) :

    x = strength

    # Change how likely it is that TSS will be applied depending on strength of core
    number = random()
    prob_sub =  {'VH' : 0.67, 'H' : 0.5, 'M' : 0.25, 'L' : 0.67}
    
    # Now write a chance of subsituting in the TSS at the expected TSS (right at the core-utr junction)
    if number <= prob_sub[x] :
        
        tssu = tss_upstream_chooser(x)
        tss_el = tss_el_chooser(x)
    
        tss = tss[:42]+tssu+tss_el
    
    return tss

def tss_upstream_chooser(strength) :

    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS upstream has three positively correlated elements and 1 negative
    tss_upstream_list = ['TTTT', 'TTCT', 'CTTA', 'AGCG']
        
    tssuD = {'VH' : [0.33, 0.67, 0.67], 'H' : [0.3, 0.6, 0.7], 'M' : [0.25, 0.5, 0.75], 'L' : [0.1, 0.2, 0.9]}

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

def tss_el_chooser(strength) :
             
    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS elements are all positively correlated, but (0) and (1) are more strongly correlated.
    # Therefore I decreased the likelihood of choosing those elements as strength decreased.
    tss_el_list = ['CAAA', 'CAAT', 'CACC', 'ACAA']
             
    tss_elD = {'VH' : [0.5, 0.5, 0.5], 'H' : [0.4, 0.5, 0.6], 'M' : [0.25, 0.5, 0.75], 'L' : [0.1, 0.5, 0.9]}

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


if __name__ == "__maine__" :
    maine()
