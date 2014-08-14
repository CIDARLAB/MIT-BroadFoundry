"""
UASGen_X
===========
    UASGen_X is the main function for creating UAS regions of promoters.
    This program reads in settings in ProGenie_Parameters.xlsx and calls SeqGen_5 for
    each UAS.  Finally, it wipes any TypeIIS sites from the sequence
    and adds on flanking sequences for cloning.    
"""

from common_functions import *
from SeqGen_5 import seqgen
from get_data import get_parameters
from site_eraser import re_eraser
from reverser import rev_comp
from uas_element_subber import *
    
def maine() :
    
    # This function pulls the data from ProGenie_Parameters.xlsx
    # I invoke above the following for loops so that the dictionary is only created once.
    parameterD = get_parameters()
    
    uasgen(parameterD)
    
def uasgen(parameterD):

    # This dictionary will hold all final sequences and metadata concerning
    # substitution
    recordD = {}
                
    # This for loop kicks off sequence generation
    for y in parameterD['uas'] :

        # Empty dictionary for sequence persistence after the loops
        uasD = {}

        # Add key to record dictionary
        recordD[y] = {}
        
        for strength in parameterD['strengths'] :
            
            # Generate random sequences with nucleotide percentages common in promoters
            seq_list = generate_uas_sequences(parameterD)

            uas_subs = uas_subber(seq_list, strength, y, parameterD)

            # UAS dictionary defined earlier, output in FASTA format,
            # and analyze the sequences with the analysis suite.          
            uasD[strength] = uas_subs[0]
            recordD[y][strength] = uas_subs[1]

        fastaD_out(uasD, 'uas%(num)s' % {'num': y})

    return recordD
    
def generate_uas_sequences(parameterD) :
    
    dataD = {'A' : parameterD['nuc_pct']['uas']['A'],
             'T' : parameterD['nuc_pct']['uas']['T'],
             'C' : parameterD['nuc_pct']['uas']['C'],
             'G' : parameterD['nuc_pct']['uas']['G'],
             'tol' : parameterD['tolerance']['uas'],
             'len' : parameterD['length']['uas'],
             'num' : parameterD['loop']['seq_number']}

    # Generate sequences for each subelement.
    # Seqgen saves theses sequences in a text file, which the next
    # lines of code will read.
    seqgen(dataD['A'], dataD['T'], dataD['C'], dataD['G'],
           dataD['tol'], dataD['len'], dataD['num'])

    # Retrieves the generated sequences from the output text file and
    # makes a list of the sequences without the names added in seqgen()
    annotated_seq_list = lister('seqsgen.txt')
            
    # Eliminate the generic FASTA names applied by seqgen
    seq_list = [line for line in annotated_seq_list if '>' not in line]

    return seq_list

def uas_subber(seq_list, strength, uas, parameterD) :

    uas_sub_recordD = {}
    
    y = uas
    
    # Following for loop substitutes the TF and poly dA:dT motifs into the promoter
    for n, seq in enumerate(seq_list) :

        uas_sub_recordD[n] = {}

        # Following logic allows some sequences to escape substitution
        if n < parameterD['loop']['no_sub_cutoff'] :
        
            tf_sub_list = tf_sub(y,
                                 strength,
                                 seq,
                                 parameterD)

            uas_sub_recordD[n]['tf_sub'] = tf_sub_list[1]
    
            pAT_sub_list = polyAT_sub(y,
                                      strength,
                                      tf_sub_list[0],
                                      parameterD)
                
            uas_sub_recordD[n]['pAT_sub'] = pAT_sub_list[1]
                
            seq_list[n] = uas_formatter(y,
                                        pAT_sub_list[0],
                                        parameterD)

            uas_sub_recordD[n]['sequence'] = seq_list[n]
        
        else :
            seq_list[n] = uas_formatter(y, seq, parameterD)
            uas_sub_recordD[n]['tf_sub'] = []
            uas_sub_recordD[n]['pAT_sub'] = []
            uas_sub_recordD[n]['sequence'] = seq_list[n]

    return [seq_list, uas_sub_recordD]

def uas_formatter(uas, seq, parameterD) :

    uas_scars = {1 : [parameterD['scars']['J2'],
                      parameterD['scars']['J3']],
                 2 : [parameterD['scars']['A'],
                      parameterD['scars']['J4']]}

    uas_flanks = {1: [parameterD['flanks']['UAS1_F'],
                      parameterD['flanks']['UAS1_R']],
                  2: [parameterD['flanks']['UAS2_F'],
                      parameterD['flanks']['UAS2_R']]}

    clean_seq = re_eraser(uas_scars[uas][0]+
                          seq+
                          uas_scars[uas][1],
                          parameterD)

    seq = uas_flanks[uas][0]+clean_seq+uas_flanks[uas][1]

    return seq
            
if __name__ == "__maine__" :
    maine()
