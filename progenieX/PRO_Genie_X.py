"""
PRO_Genie_X
===========
    PRO_Genie_X is a greedy constraint-based algorithm for generating
    nucleotide sequences that have the same features as yeast promoters.
    It is composed of a variety of functions, which are distributed among
    various python scripts according to the task that function performs.

    Running this function, PRO_Genie_X.py, will execute all other scripts
    to create a specified number of promoters according to the parameters
    set in ProGenie_Parameters.xlsx.  Therefore, the user need only set
    the desired parameters in Excel and then run this function to create
    promoters!!

    For more information, see ABOUT_PROGENIE.txt.
"""

#   Yeast Promoter Generation in Python
#   Copyright (C) 2014 by
#   Eric M. Young <ericyoung7@gmail.com>
#   All rights reserved.

__author__ = 'Eric M. Young <ericyoung7@gmail.com>, Voigt Lab, MIT'
__version__ = '2.0'

###########################################################################
# MAIN SCRIPT FOR CREATING YEAST PROMOTERS
###########################################################################

import pickle
from common_functions import *
from CoreGen_X import coregen
from UASGen_X import uasgen
from get_data import get_parameters

def maine() :
    pro_genie(int(raw_input("Number of sequences (must be multiple of 4):")))

def pro_genie(n) :

    # Get data from ProGenie_parameters.xlsx
    parameterD = get_parameters()
    
    # Generate empty record dictionary
    recordD = {}
    
    # Generate sequences using the coregen and uasgen functions
    recordD['core'] = coregen(n, parameterD)
    recordD['uas'] = uasgen(n, parameterD)

    # This is a stitching function for combining the core, uas1, and uas2
    # lists into full length promoters for analysis.  
    if False :
        promoter_stitcher(n)

    # This kicks off the regex based analysis of generated sequences
    if False :
        promoter_analysis()
    
    with open('promoter_pickle.txt', 'w') as p :
        pickle.dump(recordD, p)

    
    # This is the output function into a human-readable text file
    dictionary_output(recordD, n)

def dictionary_output(recordD, n) :

    # Define number of strengths
    n_in_strength = n/4

    strength = ['VH', 'H', 'M', 'L']
    multiplierD = {'VH': 0, 'H': 1, 'M': 2, 'L': 3}
    number = []
    
    for i in range(n_in_strength) :
        number.append(i)
    
    core_subL = ['tss_sub', 'kozak_sub', 'polyAT_sub', 'tata_sub']
    core_subD = {'tss_sub': 'TSS', 'kozak_sub': 'KOZ', 'polyAT_sub': 'pAT', 'tata_sub': 'TATAWAWR'}

    uas_subL = ['tf_sub', 'pAT_sub']    
    uasL = [1,2]

    record_filename = 'readable_record.txt'
    clear(record_filename)
    record = open(record_filename, 'a')
    
    for x in strength :
        for i in number :
            record.write('>EMYcore_%(number)s\n' % {'number' : i+1+multiplierD[x]*4})
            record.write('>>%(x)s\n' % {'x' : x})
            
            for y in core_subL :
                if len(recordD['core'][x][i][y]) > 1: 
                    record.write('>>%(name)s %(loc)s %(seq)s\n' % {'name' : core_subD[y],
                                                          'loc' : recordD['core'][x][i][y][1],
                                                          'seq' : recordD['core'][x][i][y][0]})
            record.write('%(seq)s\n' % {'seq' : recordD['core'][x][i]['sequence']})

    for u in uasL :
        for x in strength :
            for i in number :
                record.write('>EMYuas%(u)s_%(number)s\n' % {'u': u, 'number' : i+1+multiplierD[x]*4})
                record.write('>>%(x)s\n' % {'x' : x})
                
                for y in uas_subL :
                    if y is 'tf_sub' :
                        for k in recordD['uas'][u][x][i]['tf_sub'] :
                           record.write('>>%(name)s_%(site#)s %(loc1)s %(loc2)s %(site_seq)s\n' %\
                            {'name'  : recordD['uas'][u][x][i][y][k]['name'][0],
                             'site#' : recordD['uas'][u][x][i][y][k]['name'][1],
                             'loc1'  : recordD['uas'][u][x][i][y][k]['location'][0],
                             'loc2'  : recordD['uas'][u][x][i][y][k]['location'][1],
                             'site_seq' : recordD['uas'][u][x][i][y][k]['name'][2]})
                        
                    if y is 'pAT_sub' :
                        for k in recordD['uas'][u][x][i][y] :
                           record.write('>>%(name)s %(loc1)s %(loc2)s %(site_seq)s\n' %\
                            {'name'  : recordD['uas'][u][x][i][y][k]['name'][0],
                             'loc1'  : recordD['uas'][u][x][i][y][k]['location'][0],
                             'loc2'  : recordD['uas'][u][x][i][y][k]['location'][1],
                             'site_seq' : recordD['uas'][u][x][i][y][k]['name'][1]})
                            
                record.write('%(seq)s\n' % {'seq' : recordD['uas'][u][x][i]['sequence']})

    record.close()

    ''' Record dictionary structure:
    core
        strength
            number
                sequence
                tss_sub
                    [site, location]
                kozak_sub
                    [site, location]
                polyAT_sub
                    [site, location]
                tata_sub
                    [site, location]
                utr_nn_erase
                tss_nn_erase
                atg_erase

    UAS
        UAS number
            strength
                number
                    sequence
                    tf_sub
                        number
                            name
                                [name, sequence variant number, sequence]
                            location
                                [start, stop]
                    pAT_sub
                        number
                            name
                                [name, sequence]
                            location
                                [start, stop]                         
    '''

def promoter_stitcher(n) :

    # Retrieve generated sequences from the output txt files.  Omit TypeIIS sites on ends
    # to enable stitching as they would be in a TypeIIS reaction.
    corelist = lister('coregen.txt')

    cores = [line[11:145] for line in corelist if '>' not in line]
    x = [line[9:10] for line in corelist if '>' in line]

    uas1list = lister('uas1gen.txt')

    uas1s = [line[11:141] for line in uas1list if '>' not in line]
    
    uas2list = lister('uas2gen.txt')

    uas2s = [line[11:141] for line in uas2list if '>' not in line]

    # Stitch together the list of sequences into promoters. Note that this is not random stitching,
    # rather sequences are combined in the order they are generated.
    pros = cores

    for n, pro in enumerate(pros) :
        pro = uas2s[n]+uas1s[n]+cores[n]
        pros[n] = pro

    # Name output data files, format sequence dictionary into FASTA, and analyze for motifs.
    prefix = 'pro'
    fasta_out(pros, prefix)
    
if __name__ == "__maine__" :
    maine()
