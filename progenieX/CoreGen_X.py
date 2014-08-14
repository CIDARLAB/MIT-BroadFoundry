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
from site_eraser import *
from core_element_subber import *

def maine() :
    
    # This function pulls the data from ProGenie_Parameters.xlsx
    # I invoke above the following for loops so that the dictionary is only created once.
    parameterD = get_parameters()
    
    coregen(parameterD)

def coregen(parameterD) :
    
    # This dictionary is an empty dictionary into which the generated sequences will
    # be substituted
    subelementD = {'VH': {'tbp': [], 'tss': [], 'utr': []},
                   'H': {'tbp': [], 'tss': [], 'utr': []},
                   'M': {'tbp': [], 'tss': [], 'utr': []},
                   'L': {'tbp': [], 'tss': [], 'utr': []}}
    
    # This dictionary will hold all final sequences and metadata concerning substitution
    recordD = {}

    # This for loop kicks of sequence generation
    for strength in parameterD['strengths'] :
        
        for subpart in parameterD['subpart'] :

            # Call the sequence generation function that handles calling SeqGen_5
            seq_list = generate_core_sequences(strength, subpart, parameterD)

            # Begin substitution logic.  Calls function for each subelement
            if subpart is 'tbp':

                tbp_subs = tbp_subber(seq_list, strength, parameterD)

                seq_list, tbp_sub_recordD = tbp_subs[0], tbp_subs[1]

            if subpart is 'tss' :

                tss_subs = tss_subber(seq_list, strength, parameterD)

                seq_list, tss_sub_recordD = tss_subs[0], tss_subs[1]
                
            if subpart is 'utr':

                utr_subs = utr_subber(seq_list, strength, parameterD)

                seq_list, utr_sub_recordD = utr_subs[0], utr_subs[1]
                    
            # This populates the different segments into the empty dictionary to hold
            # all sequences together.
            subelementD[strength][subpart] = seq_list

        # Combine all record dictionaries into one record dictionary
        recordD[strength] = {}
        
        for n, seq in enumerate(seq_list) :
                
            recordD[strength][n] = {}
                
            for k in tbp_sub_recordD[n] :

                recordD[strength][n][k] = tbp_sub_recordD[n][k]

            for k in tss_sub_recordD[n] :

                recordD[strength][n][k] = tss_sub_recordD[n][k]

            for k in utr_sub_recordD[n] :

                recordD[strength][n][k] = utr_sub_recordD[n][k]

    # Stitch subelements
    coreD = stitch_subelements(parameterD, subelementD)
    
    # Format sequence dictionary into FASTA and save to text file. 
    fastaD_out(coreD, 'core')

    # Add final sequences to the record dictionary.
    for strength in parameterD['strengths'] :
        for n in range(parameterD['loop']['seq_number']) :
            recordD[strength][n]['sequence'] = coreD[strength][n]
    
    return recordD

def generate_core_sequences(strength, subpart, parameterD) :
    
    dataD = {'A' : parameterD['nuc_pct'][strength][subpart]['A'],
             'T' : parameterD['nuc_pct'][strength][subpart]['T'],
             'C' : parameterD['nuc_pct'][strength][subpart]['C'],
             'G' : parameterD['nuc_pct'][strength][subpart]['G'],
             'tol' : parameterD['tolerance'][subpart],
             'len' : parameterD['length'][subpart],
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

def tbp_subber(seq_list, strength, parameterD):

    tbp_sub_recordD = {}
    
    # This loop substitutes TATA and polyAT sites into the tbp subelement.
    # This arrangement allows nucleosome free region close to the tss.
    
    for n, seq in enumerate(seq_list):

        tbp_sub_recordD[n] = {}

        # Following logic allows some sequences to escape substitution
        if n < parameterD['loop']['no_sub_cutoff'] :
                        
            tata_sub_list = tata_sub(seq,
                                     strength,
                                     parameterD)

            tbp_sub_recordD[n]['tata_sub'] = tata_sub_list[1]

            polyAT_sub_list = polyAT_tbp_sub(tata_sub_list[0],
                                             strength,
                                             parameterD)

            tbp_sub_recordD[n]['polyAT_sub'] = polyAT_sub_list[1]
                        
            seq_list[n] = polyAT_sub_list[0]
                        
        else :
            tbp_sub_recordD[n]['tata_sub'] = []
            tbp_sub_recordD[n]['polyAT_sub'] = []


    return [seq_list, tbp_sub_recordD]

def tss_subber(seq_list, strength, parameterD) :
    
    # Following block substitutes transcription start sites from Lubliner
    # et al. as well as erases NAB3 and NRD1 mRNA degradation signals in
    # the forward direction.  This is to reduce likelihood of creating
    # an easily degraded 5'-UTR.

    tss_sub_recordD = {}
    
    for n, seq in enumerate(seq_list) :

        tss_sub_recordD[n] = {}
        
        # Following logic allows some sequences to escape substitution
        if n < parameterD['loop']['no_sub_cutoff'] :
                        
            tss_sub_list = tss_sub(seq,
                                   strength,
                                   parameterD)
                    
            tss_sub_recordD[n]['tss_sub'] = tss_sub_list[1]
                    
            tss_nn_erase_list = nab_nrd_eraser(tss_sub_list[0],
                                               strength,
                                               parameterD)

            tss_sub_recordD[n]['tss_nn_erase'] = tss_nn_erase_list[1]

            seq_list[n] = tss_nn_erase_list[0]
                        
        else : 
            tss_sub_recordD[n]['tss_sub'] = []
            tss_sub_recordD[n]['tss_nn_erase'] = []

    return [seq_list, tss_sub_recordD]
                        
def utr_subber(seq_list, strength, parameterD) :

    # Following block substitutes consensus kozak sequences from Dvir et al.
    # as well as erases NAB3 and NRD1 mRNA degradation signals in the forward
    # direction.  This is to reduce likelihood of creating an easily degraded
    # 5'-UTR. Furthermore, since upstream ATG can reduce expression, also
    # erase any ATG that may have arisen.

    utr_sub_recordD = {}
    
    for n, seq in enumerate(seq_list) :

        utr_sub_recordD[n] = {}
        
        # Following logic allows some sequences to escape substitution
        if n < parameterD['loop']['no_sub_cutoff'] :
                        
            kozak_sub_list = kozak_sub(seq,
                                       strength,
                                       parameterD)
                    
            utr_sub_recordD[n]['kozak_sub'] = kozak_sub_list[1]
                    
            utr_nn_erase_list = nab_nrd_eraser(kozak_sub_list[0],
                                               strength,
                                               parameterD)

            utr_sub_recordD[n]['utr_nn_erase'] = utr_nn_erase_list[1]

            atg_erase_list = atg_eraser(utr_nn_erase_list[0],
                                        strength,
                                        parameterD)

            utr_sub_recordD[n]['atg_erase'] = atg_erase_list[1]
                    
            seq_list[n] = atg_erase_list[0]

        else :
            utr_sub_recordD[n]['kozak_sub'] = []
            utr_sub_recordD[n]['utr_nn_erase'] = []
            utr_sub_recordD[n]['atg_erase'] =  []

    return [seq_list, utr_sub_recordD]

def stitch_subelements(parameterD, subelementD) :
    
    # Define empty dictionary to substitute stitched sequences
    coreD = {'VH': [], 'H': [], 'M': [], 'L': []}

    # This stitches together the subsegments into a full core promoter and adds TypeIIS sites.
    for strength in parameterD['strengths'] :
        
        coreD[strength] = subelementD[strength]['tbp']
        
        for n, seq in enumerate(subelementD[strength]['tbp']):
            
            core = re_eraser(parameterD['scars']['J1']+
                             subelementD[strength]['tbp'][n]+
                             subelementD[strength]['tss'][n]+
                             subelementD[strength]['utr'][n]+
                             parameterD['scars']['B'],
                             parameterD)

            coreD[strength][n] = parameterD['flanks']['Core_F']+core+parameterD['flanks']['Core_R']

    return coreD

if __name__ == "__maine__" :
    maine()
