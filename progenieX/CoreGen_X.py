"""
CoreGen_X
===========
    CoreGen_X is the main function for creating segments of core promoters.
    This program reads in settings in ProGenie_Parameters.xlsx and calls SeqGen_5 for
    each subelement of the core promoter: tata-binding region (tbp), transcription start
    site (tss), and anticipated 5'-UTR (utr).  It then passes these sequences to the
    various substitution functions.  Finally, it wipes any TypeIIS sites from the sequence
    and adds on flanking sequences for cloning.    
"""

from common_functions import *
from SeqGen_5 import seqgen
from get_data import get_parameters
from core_element_subber import *
from site_eraser import *

def maine() :
    
    # This function pulls the data from ProGenie_Parameters.xlsx
    # I invoke above the following for loops so that the dictionary is only created once.
    parameterD = get_parameters()
    
    coregen(int(raw_input("Number (must be a multiple of 4):")), parameterD)

def coregen(core_number, parameterD) :

    # Common iterating list for sequence scripts
    ATCG = ['A', 'T', 'C', 'G']
    
    # Define iterating lists for strengths and each core subpart 
    strengths = ['VH', 'H', 'M', 'L']
    subpart = ['tbp', 'tss', 'utr']

    # This dictionary is an empty dictionary into which the generated sequences will
    # be substituted
    subelementD = {'VH': {'tbp': [], 'tss': [], 'utr': []},
                   'H': {'tbp': [], 'tss': [], 'utr': []},
                   'M': {'tbp': [], 'tss': [], 'utr': []},
                   'L': {'tbp': [], 'tss': [], 'utr': []}}

    # Since the generator will create VH, H, L and M for each cornum, to output the
    # number of desired sequences, the input has to be divided by 4.
    cores = core_number/4

    # This dictionary will hold all final sequences and metadata concerning substitution
    recordD = {}
    for x in strengths :
        recordD[x] = {}
        for n in range(cores) :
            recordD[x][n] = {}

    # This for loop kicks of sequence generation
    for x in strengths :
        
        for y in subpart :

            # Call the sequence generation function that handles calling SeqGen_5 using the
            # parameters set in ProGenie_Parameters.xlsx
            seq_list = generate_core_sequences(cores, x, y, parameterD)
         
            # This loop substitutes TATA and polyAT sites into the tbp region.
            # This arrangement allows for one more pdW, but without combinatorial
            # explosion in UAS1, and also allows nucleosome free region close to the
            # tss. 
            if y is 'tbp':
                for n, seq in enumerate(seq_list):
                        
                    tata_sub_list = tata_sub(seq, x)

                    recordD[x][n]['tata_sub'] = tata_sub_list[1]

                    polyAT_sub_list = polyAT_tbp_sub(tata_sub_list[0],x)

                    recordD[x][n]['polyAT_sub'] = polyAT_sub_list[1]
                        
                    seq_list[n] = polyAT_sub_list[0]
                    
            # Since there is data available about elements around the TSS from Lubliner
            # I also wrote a function that puts in TSS elements right at the end of the
            # tss element
            if y is 'tss' :
                for n, seq in enumerate(seq_list) :
                    
                    tss_sub_list = tss_sub(seq, x)
                    
                    recordD[x][n]['tss_sub'] = tss_sub_list[1]
                    
                    tss_nn_erase_list = nab_nrd_eraser(tss_sub_list[0], x)

                    recordD[x][n]['tss_nn_erase'] = tss_nn_erase_list[1]

                    seq_list[n] = tss_nn_erase_list[0]
                    
            # Since the consensus Kozak is also unlikely to appear at the end of the UTR,
            # I also wrote a function that substitutes it in.
            if y is 'utr':
                for n, seq in enumerate(seq_list) :
                    
                    kozak_sub_list = kozak_sub(seq, x)
                    
                    recordD[x][n]['kozak_sub'] = kozak_sub_list[1]
                    
                    atg_erase_list = atg_eraser(kozak_sub_list[0], x)

                    recordD[x][n]['atg_erase'] = atg_erase_list[1]
                    
                    utr_nn_erase_list = nab_nrd_eraser(atg_erase_list[0], x)

                    recordD[x][n]['utr_nn_erase'] = utr_nn_erase_list[1]
                    
                    seq_list[n] = utr_nn_erase_list[0]
                    
            # This populates the different segments into the empty dictionary to hold
            # all sequences together.
            subelementD[x][y] = seq_list

    # Stitch subelements
    coreD = stitch_subelements(parameterD, subelementD, strengths)
    
    # Format sequence dictionary into FASTA and save to text file. 
    fastaD_out(coreD, 'core')

    # Add final sequences to the record dictionary.
    for x in strengths :
        for n in range(cores) :
            recordD[x][n]['sequence'] = coreD[x][n]
    
    return recordD

def generate_core_sequences(cores, x, y, parameterD) :
    
    dataD = {'A' : parameterD['nuc_pct'][x][y]['A'],
             'T' : parameterD['nuc_pct'][x][y]['T'],
             'C' : parameterD['nuc_pct'][x][y]['C'],
             'G' : parameterD['nuc_pct'][x][y]['G'],
             'tol' : parameterD['tolerance'][y],
             'len' : parameterD['length'][y],
             'num' : cores}

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

def stitch_subelements(parameterD, subelementD, strengths) :
    
    # Define empty dictionary to substitute stitched sequences
    coreD = {'VH': [], 'H': [], 'M': [], 'L': []}

    # This stitches together the subsegments into a full core promoter and adds TypeIIS sites.
    for x in strengths :
        
        coreD[x] = subelementD[x]['tbp']
        
        for n, seq in enumerate(subelementD[x]['tbp']):
            
            core = re_eraser(parameterD['scars']['J1']+
                             subelementD[x]['tbp'][n]+
                             subelementD[x]['tss'][n]+
                             subelementD[x]['utr'][n]+
                             parameterD['scars']['B'])

            coreD[x][n] = parameterD['flanks']['Core_F']+core+parameterD['flanks']['Core_R']

    return coreD

if __name__ == "__maine__" :
    maine()
