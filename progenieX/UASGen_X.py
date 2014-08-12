from xlrd import *
from common_functions import *
from excel_functions import cell
from reverser import rev_comp
from SeqGen_5 import seqgen
from get_data import get_parameters
from site_eraser import re_eraser
from uas_subs import *
    
def maine() :
    
    # This function pulls the data from ProGenie_Parameters.xlsx
    # I invoke above the following for loops so that the dictionary is only created once.
    parameterD = get_parameters()
    
    uasgen(int(raw_input("Number (must be a multiple of 4):")), parameterD)
    
def uasgen(uas_number, parameterD):

    # Define names of strength levels and iterable lists
    strengths = ['VH', 'H', 'M', 'L']
    ATCG = ['A', 'T', 'C', 'G']
    uas_num = [1,2]

    uas_scars = {1 : [parameterD['scars']['J2'],
                      parameterD['scars']['J3']],
                 2 : [parameterD['scars']['A'],
                      parameterD['scars']['J4']]}

    uas_flanks = {1: [parameterD['flanks']['UAS1_F'],
                      parameterD['flanks']['UAS1_R']],
                  2: [parameterD['flanks']['UAS2_F'],
                      parameterD['flanks']['UAS2_R']]}

    # Since the generator will create VH, H, L and M for each uasnum, to output the
    # number of desired sequences, the input has to be divided by 4.
    uas_n = uas_number/4
    
    # This dictionary will hold all final sequences and metadata concerning substitution
    recordD = {}
    for y in uas_num :
        recordD[y] = {}
        for x in strengths :
            recordD[y][x] = {}
            for n in range(uas_n) :
                recordD[y][x][n] = {}
                
    # This for loop kicks off sequence generation
    for y in uas_num :
        
        uasD = {}
        
        for x in strengths :
            
            # Generate random sequences with nucleotide percentages common in promoters
            seq_list = generate_uas_sequences(uas_n, parameterD)

            # This for loop substitutes the TF and poly dA:dT motifs into the promoter
            for n, seq in enumerate(seq_list) :
                
                tf_sub_list = tf_sub(y, x, seq, parameterD)

                recordD[y][x][n]['tf_sub'] = tf_sub_list[1]

                pAT_sub_list = polyAT_sub(y, x, tf_sub_list[0], parameterD)
                
                recordD[y][x][n]['pAT_sub'] = pAT_sub_list[1]
                
                uas = re_eraser(uas_scars[y][0]+pAT_sub_list[0]+uas_scars[y][1])
       
                seq_list[n] = uas_flanks[y][0]+uas+uas_flanks[y][1]

                recordD[y][x][n]['sequence'] = seq_list[n]

            # UAS dictionary defined earlier, output in FASTA format,
            # and analyze the sequences with the analysis suite.          
            uasD[x] = seq_list        

        fastaD_out(uasD, 'uas%(num)s' % {'num': y})

    return recordD
    
def generate_uas_sequences(uas_n, parameterD) :
    
    dataD = {'A' : parameterD['nuc_pct']['uas']['A'],
             'T' : parameterD['nuc_pct']['uas']['T'],
             'C' : parameterD['nuc_pct']['uas']['C'],
             'G' : parameterD['nuc_pct']['uas']['G'],
             'tol' : parameterD['tolerance']['uas'],
             'len' : parameterD['length']['uas'],
             'num' : uas_n}

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
            
if __name__ == "__maine__" :
    maine()
