import re
from itertools import product
from xlrd import *
from reverser import rev_comp

def maine():
    get_parameters()
    
def get_parameters():
    
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'

    # Define number of sequences and number of unsubstituted sequences
    # Desired
    inputD = {'total'  : int(cell(iB,gS,'C1')),
              'no_sub' : int(cell(iB,gS,'C2'))}

    # Define numbers to be used in loops, since loops go by four strength
    # categories, the total number needs to be divided by four.
    loopD = {'seq_number'    :  inputD['total']/4,
             'no_sub_cutoff' : (inputD['total']-inputD['no_sub'])/4}
    
    #################
    #Iterating Lists#
    #################

    # Common iterating list for sequence scripts
    ATCG = ['A', 'T', 'C', 'G']
    
    # Define iterating lists for strengths and each core subpart & UAS 
    strengths = ['VH', 'H', 'M', 'L']
    subpart = ['tbp', 'tss', 'utr']
    uas_num = [1,2]


    ################################
    #Sequence Generation Parameters#
    ################################

    # Associated fractions from Lubliner et al. 2013 for the core elements
    # ("tbp", "tss", and "utr") and the uas elements
    nuc_pctD = {'VH': {'tbp': {'A': cell(iB,cS,'C3'), 'T': cell(iB,cS,'D3'),
                               'C': cell(iB,cS,'E3'), 'G': cell(iB,cS,'F3')},
                       'tss': {'A': cell(iB,cS,'C4'), 'T': cell(iB,cS,'D4'),
                               'C': cell(iB,cS,'E4'), 'G': cell(iB,cS,'F4')},
                       'utr': {'A': cell(iB,cS,'C5'), 'T': cell(iB,cS,'D5'),
                               'C': cell(iB,cS,'E5'), 'G': cell(iB,cS,'F5')}},
                'H' : {'tbp': {'A': cell(iB,cS,'C6'), 'T': cell(iB,cS,'D6'),
                              'C': cell(iB,cS,'E6'), 'G': cell(iB,cS,'F6')},
                       'tss': {'A': cell(iB,cS,'C7'), 'T': cell(iB,cS,'D7'),
                               'C': cell(iB,cS,'E7'), 'G': cell(iB,cS,'F7')},
                       'utr': {'A': cell(iB,cS,'C8'), 'T': cell(iB,cS,'D8'),
                               'C': cell(iB,cS,'E8'), 'G': cell(iB,cS,'F8')}},
                'M' : {'tbp': {'A': cell(iB,cS,'C9'), 'T': cell(iB,cS,'D9'),
                               'C': cell(iB,cS,'E9'), 'G': cell(iB,cS,'F9')},
                       'tss': {'A': cell(iB,cS,'C10'), 'T': cell(iB,cS,'D10'),
                               'C': cell(iB,cS,'E10'), 'G': cell(iB,cS,'F10')},
                       'utr': {'A': cell(iB,cS,'C11'), 'T': cell(iB,cS,'D11'),
                               'C': cell(iB,cS,'E11'), 'G': cell(iB,cS,'F11')}},
                'L' : {'tbp': {'A': cell(iB,cS,'C12'), 'T': cell(iB,cS,'D12'),
                              'C': cell(iB,cS,'E12'), 'G': cell(iB,cS,'F12')},
                       'tss': {'A': cell(iB,cS,'C13'), 'T': cell(iB,cS,'D13'),
                               'C': cell(iB,cS,'E13'), 'G': cell(iB,cS,'F13')},
                       'utr': {'A': cell(iB,cS,'C14'), 'T': cell(iB,cS,'D14'),
                               'C': cell(iB,cS,'E14'), 'G': cell(iB,cS,'F14')}},
                'uas' : {'A' : int(cell(iB,uS,'A3')), 'T' : int(cell(iB,uS,'B3')),
                         'C' : int(cell(iB,uS,'C3')), 'G' : int(cell(iB,uS,'D3'))}}

    # Tolerance for variation around above set nucleotide percentages
    toleranceD = {'tbp': cell(iB,cS,'G3'),
                  'tss': cell(iB,cS,'G4'),
                  'utr': cell(iB,cS,'G5'),
                  'uas' : float(cell(iB,uS,'E3'))}

    # Length of each subelement
    lengthD = {'tbp': int(cell(iB,cS,'H3')),
               'tss': int(cell(iB,cS,'H4')),
               'utr': int(cell(iB,cS,'H5')),
               'uas': int(cell(iB,uS,'F3'))}

    ######################################
    #Core Element Substitution Parameters#
    ######################################

    # How many sites have a chance to add for core elements
    core_siteD = {'pAT' : int(cell(iB,cS,'B43'))}
    
    # Slice Locations - called "Insertion Locations" in ProGenie_Parameters
    core_sliceD = {'tata'   : [int(cell(iB,cS,'C18')),
                               int(cell(iB,cS,'C19')),
                               int(cell(iB,cS,'C20'))],
                   'tss'    :  int(cell(iB,cS,'C21')),
                   'kozak'  : [int(cell(iB,cS,'C22')),
                               int(cell(iB,cS,'C23'))],
                   'tbp_pAT':  int(cell(iB,cS,'C24'))}

    # Probability of substitution dictionary for core elements
    core_sub_probD = {'tata' : [cell(iB, cS, 'B26'),
                                cell(iB, cS, 'D28'),
                                cell(iB, cS, 'D29')],
                      'tss'  : {'VH': cell(iB, cS, 'C36'),
                                'H' : cell(iB, cS, 'C37'),
                                'M' : cell(iB, cS, 'C38'),
                                'L' : cell(iB, cS, 'C39')},
                      'kozak': {'VH': cell(iB, cS, 'C31'),
                                'H' : cell(iB, cS, 'C32'),
                                'M' : cell(iB, cS, 'C33'),
                                'L' : cell(iB, cS, 'C34'),
                                'db': cell(iB, cS, 'C35')},
                      'tbp_pAT' : {'VH': float(cell(iB,cS,'B46')),
                                   'H' : float(cell(iB,cS,'C46')),
                                   'M' : float(cell(iB,cS,'D46')),
                                   'L' : float(cell(iB,cS,'E46'))}}
    
    # Kozak mutants derived from consensus sequences in Dvir et al. In the design,
    # the -1 nucleotide will always be A because Bscar is AATG.
    # So the Kozaks are only 9 nucleotides long. The first two Kozaks
    # are positively correlated, the second two are negatively correlated.
        
    kozak_seqL = [cell(iB, cS, 'R19'),
                  cell(iB, cS, 'R20'),
                  cell(iB, cS, 'R21'),
                  cell(iB, cS, 'R22')]
        
    kozak_choiceD = {'VH': [cell(iB, cS, 'I20'),
                            cell(iB, cS, 'I21'),
                            cell(iB, cS, 'I22')],
                     'H' : [cell(iB, cS, 'K20'),
                            cell(iB, cS, 'K21'),
                            cell(iB, cS, 'K22')],
                     'M' : [cell(iB, cS, 'M20'),
                            cell(iB, cS, 'M21'),
                            cell(iB, cS, 'M22')],
                     'L' : [cell(iB, cS, 'O20'),
                            cell(iB, cS, 'O21'),
                            cell(iB, cS, 'O22')]}
    
    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS upstream has three positively correlated elements and 1 negative
    tss_upstreamL = [cell(iB,cS,'R28'),
                     cell(iB,cS,'R29'),
                     cell(iB,cS,'R30'),
                     cell(iB,cS,'R31')]
      
    tss_upstream_choiceD = {'VH' : [cell(iB,cS,'I29'),
                                    cell(iB,cS,'I30'),
                                    cell(iB,cS,'I31')],
                            'H'  : [cell(iB,cS,'K29'),
                                    cell(iB,cS,'K30'),
                                    cell(iB,cS,'K31')],
                            'M'  : [cell(iB,cS,'M29'),
                                    cell(iB,cS,'M30'),
                                    cell(iB,cS,'M31')],
                            'L'  : [cell(iB,cS,'O29'),
                                    cell(iB,cS,'O30'),
                                    cell(iB,cS,'O31')]}

    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS elements are all positively correlated, but (0) and (1) are more strongly correlated.
    # Therefore I decreased the likelihood of choosing those elements as strength decreased.
    tss_elL = [cell(iB,cS,'R33'),
               cell(iB,cS,'R34'),
               cell(iB,cS,'R35'),
               cell(iB,cS,'R36')]
             
    tss_el_choiceD = {'VH': [cell(iB,cS,'I34'),
                             cell(iB,cS,'I35'),
                             cell(iB,cS,'I36')],
                     'H'  : [cell(iB,cS,'K34'),
                             cell(iB,cS,'K35'),
                             cell(iB,cS,'K36')],
                     'M'  : [cell(iB,cS,'M34'),
                             cell(iB,cS,'M35'),
                             cell(iB,cS,'M36')],
                     'L'  : [cell(iB,cS,'O34'),
                             cell(iB,cS,'O35'),
                             cell(iB,cS,'O36')]}
                    
    # Probability of choosing a pAT sequence for the 'tbp' subelement
    # of the core.
    tbp_pAT_choice_probD = {'VH' : [float(cell(iB,cS,'H40')),
                                    float(cell(iB,cS,'H41'))],
                            'H'  : [float(cell(iB,cS,'J40')),
                                    float(cell(iB,cS,'J41'))],
                            'M'  : [float(cell(iB,cS,'L40')),
                                    float(cell(iB,cS,'L41'))],
                            'L'  : [float(cell(iB,cS,'N40')),
                                    float(cell(iB,cS,'N41'))]}

    # Define the potential sequences to substitute. Poly dA:dT - 1 T, 1 A, 1 random mix
    tbp_pAT_seqD = {0: ['pT' , cell(iB,cS,'R42')],
                    1: ['pA' , cell(iB,cS,'R41')],
                    2: ['mix', cell(iB,cS,'R40')]}


    #####################################
    #UAS Element Substitution Parameters#
    #####################################
                    
    # How many TFBS & polyAT to have a chance to add
    tf_siteD = {1 : int(cell(iB,uS,'A8')),
                2 : int(cell(iB,uS,'C8'))}

    polyAT_siteD = {1 : int(cell(iB,uS,'B8')),
                    2 : int(cell(iB,uS,'D8'))}

    # Where to add elements (polyAT before site & TFBS after site for UAS2
    # (polyAT after site & TFBS before site for UAS1)
    tf_lociD = {1 : [int(cell(iB,uS,'A13')), int(cell(iB,uS,'A14'))],
                2 : [int(cell(iB,uS,'C13')), int(cell(iB,uS,'C14')),
                     int(cell(iB,uS,'C15')), int(cell(iB,uS,'C16'))]}

    polyAT_lociD = {1 : [int(cell(iB,uS,'B13')), int(cell(iB,uS,'B14')),
                         int(cell(iB,uS,'B15')), int(cell(iB,uS,'B16'))],
                    2 : [int(cell(iB,uS,'D13'))]}

    # How likley TFBS addition occurs at subsequent sites past the first site
    tf_addD = {1 : {'VH': float(cell(iB,uS,'B20')), 'H' : float(cell(iB,uS,'C20')),
                    'M' : float(cell(iB,uS,'D20')), 'L' : float(cell(iB,uS,'E20'))},
               2 : {'VH': float(cell(iB,uS,'B21')), 'H' : float(cell(iB,uS,'C21')),
                    'M' : float(cell(iB,uS,'D21')), 'L' : float(cell(iB,uS,'E21'))}}

    tf_choice_probD = {1: [float(cell(iB,uS,'C25')), float(cell(iB,uS,'C26')),
                           float(cell(iB,uS,'C27')), float(cell(iB,uS,'C28')),
                           float(cell(iB,uS,'C29'))],
                       2: [float(cell(iB,uS,'E25')), float(cell(iB,uS,'E26')),
                           float(cell(iB,uS,'E27')), 1, 1]}
            
    # How likely poly dA:dT addition occurs
    pAT_addD = {1 : {'VH': float(cell(iB,uS,'B34')), 'H' : float(cell(iB,uS,'C34')),
                     'M' : float(cell(iB,uS,'D34')), 'L' : float(cell(iB,uS,'E34'))},
                2 : {'VH': float(cell(iB,uS,'B35')), 'H' : float(cell(iB,uS,'C35')),
                     'M' : float(cell(iB,uS,'D35')) ,'L' : float(cell(iB,uS,'E35'))}}
    
    pAT_choice_probD = {'VH' : [float(cell(iB,uS,'J45')), float(cell(iB,uS,'J46'))],
                        'H'  : [float(cell(iB,uS,'L45')), float(cell(iB,uS,'L46'))],
                        'M'  : [float(cell(iB,uS,'N45')), float(cell(iB,uS,'N46'))],
                        'L'  : [float(cell(iB,uS,'P45')), float(cell(iB,uS,'P46'))]}      

    # Define the potential sequences to substitute. Poly dA:dT - 1 T, 1 A, 1 random mix
    pAT_seqD = {0: ['pT',  cell(iB,uS,'S47')],
                1: ['pA',  cell(iB,uS,'S45')],
                2: ['mix', cell(iB,uS,'S46')]}

    ###########
    #TFBS Selection Parameters
    ##########

    # REB1
    # BINDING SITES FROM YEASTRACT AND MOGNO ET AL.
    reb1 = [cell(iB,uS,'S5'), cell(iB,uS,'S6'), cell(iB,uS,'S7'), cell(iB,uS,'S8')]

    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # for VH strength, only the first or second REB1 site is possible to choose.  For H strength, three are possible, adding
    # in the consensus site published in YEASTRACT.  For M, the 'weak' site from Mogno is included as a choice. For L, the
    # weak site is most likely.
    reb1_site_choiceD = {'VH' : [cell(iB,uS,'J6'), cell(iB,uS,'J7'), cell(iB,uS,'J8')],
                         'H' : [cell(iB,uS,'L6'), cell(iB,uS,'L7'), cell(iB,uS,'L8')],
                         'M' : [cell(iB,uS,'N6'), cell(iB,uS,'N7'), cell(iB,uS,'N8')],
                         'L' : [cell(iB,uS,'P6'), cell(iB,uS,'P7'), cell(iB,uS,'P8')]}

    # RAP1
    # BINDING SITES FROM YEASTRACT AND MOGNO ET AL.
    rap1 = [cell(iB,uS,'S14'), cell(iB,uS,'S15'), cell(iB,uS,'S16')]

    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # Very similar to the code for the REB1_chooser, so check that out for more discussion.
    rap1_site_choiceD = {'VH' : [cell(iB,uS,'J15'), cell(iB,uS,'J16')],
                         'H' : [cell(iB,uS,'L15'), cell(iB,uS,'L16')],
                         'M' : [cell(iB,uS,'N15'), cell(iB,uS,'N16')],
                         'L' : [cell(iB,uS,'P15'), cell(iB,uS,'P16')]}
    
    # GCR1
    # BINDING SITES FROM YEASTRACT, MOGNO ET AL., AND BAI ET AL.
    gcr1 = [cell(iB,uS,'S22'), cell(iB,uS,'S23'), cell(iB,uS,'S24'), cell(iB,uS,'S25')]
    
    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # Very similar to the code for the REB1_chooser, so check that out for more discussion.
    gcr1_site_choiceD = {'VH' : [cell(iB,uS,'J23'), cell(iB,uS,'J24'), cell(iB,uS,'J25')],
                         'H' : [cell(iB,uS,'L23'), cell(iB,uS,'L24'), cell(iB,uS,'L25')],
                         'M' : [cell(iB,uS,'N23'), cell(iB,uS,'N24'), cell(iB,uS,'N25')],
                         'L' : [cell(iB,uS,'P23'), cell(iB,uS,'P24'), cell(iB,uS,'P25')]}

    # ABF1
    # BINDING SITES FROM YEASTRACT, LAST TWO BASED OFF "RTCRYYYNNNACG" 
    abf1 = [cell(iB,uS,'S30'), cell(iB,uS,'S31'), cell(iB,uS,'S32')]

    abf1_site_choiceL = [cell(iB,uS,'J31'), cell(iB,uS,'J32')]

    # MCM1
    # BINDING SITES FROM BAI ET AL.
    mcm1 = [cell(iB,uS,'S37'), cell(iB,uS,'S38')]

    mcm1_site_choiceL = [cell(iB,uS,'J38')]

    # RSC3
    # BINDING SITE FROM BAI ET AL.
    rsc3 = cell(iB,uS,'S41')

    #########################
    #Site Removal Parameters#
    #########################

    reD = {"bpiI_f": {'seq' : cell(iB,gS,'B20'),
                      'fix' : cell(iB,gS,'D20')},
           "bpiI_r": {'seq' : rev_comp(cell(iB,gS,'B20')),
                      'fix' : cell(iB,gS,'D21')},
           "bsaI_f": {'seq' : cell(iB,gS,'B21'),
                      'fix' : cell(iB,gS,'D22')},
           "bsaI_r": {'seq' : rev_comp(cell(iB,gS,'B21')),
                      'fix' : cell(iB,gS,'D23')},
           "sapI_f": {'seq' : cell(iB,gS,'B22'),
                      'fix' : cell(iB,gS,'D24')},
           "sapI_r": {'seq' : rev_comp(cell(iB,gS,'B22')),
                      'fix' : cell(iB,gS,'D25')},
           "mlyI_f": {'seq' : cell(iB,gS,'B23'),
                      'fix' : cell(iB,gS,'D26')},
           "mlyI_r": {'seq' : rev_comp(cell(iB,gS,'B23')),
                      'fix' : cell(iB,gS,'D27')}}

    atgD = {'atg' :  cell(iB,gS,'B30'),
            'fix' : [cell(iB,gS,'B31'),
                     cell(iB,gS,'B32'),
                     cell(iB,gS,'B33')],
            'fix_choice' : [cell(iB,gS,'D32'),
                            cell(iB,gS,'D33')]}

    nnD = {'nab3'  :  cell(iB,gS,'B37'),
           'nrd1_1':  cell(iB,gS,'B38'),
           'nrd1_2':  cell(iB,gS,'B39'),
           'fix'   : [cell(iB,gS,'C37'),
                      cell(iB,gS,'C38')]}
    

    #############################
    #Cloning Sequence Parameters#
    #############################
    
    # Scars for cloning
    scarD = {'A'  : cell(iB,gS,'B5'),   # GTGC - UAS2-F
             'B'  : cell(iB,gS, 'B6'),  # AATG - CORE-R
             'J1' : cell(iB,gS,'B7'),   # TTCT - CORE-F
             'J2' : cell(iB,gS,'B8'),   # AAAC - UAS1-F
             'J3' : cell(iB,gS,'B9'),   # GACC - UAS1-R
             'J4' : cell(iB,gS,'B10')}  # CCGA - UAS2-R

    # Flanking sequences
    flankD = {'Core_F' : cell(iB,gS,'G5'),           # BsaI with primer anneal overhang
              'Core_R' : rev_comp(cell(iB,gS,'G6')), # Reverse BsaI with overhang
              'UAS1_F' : cell(iB,gS,'G7'),
              'UAS1_R' : rev_comp(cell(iB,gS,'G8')),
              'UAS2_F' : cell(iB,gS,'G9'),
              'UAS2_R' : rev_comp(cell(iB,gS,'G10'))}
    
    
    # Combine all parameters into one dictionary
    parameterD = {'ATCG' : ATCG,
                  'strengths' : strengths,
                  'subpart' : subpart,
                  'uas' : uas_num,
                  'input' : inputD,
                  'loop' : loopD,
                  'nuc_pct': nuc_pctD,
                  'tolerance' : toleranceD,
                  'length' : lengthD,
                  'core_site' : core_siteD,
                  'core_slice' : core_sliceD,
                  'core_sub_prob' : core_sub_probD,
                  'kozak' : kozak_seqL,
                  'kozak_choice' : kozak_choiceD,
                  'tss_upstream' : tss_upstreamL,
                  'tss_upstream_choice' : tss_upstream_choiceD,
                  'tss_el' : tss_elL,
                  'tss_el_choice' : tss_el_choiceD,
                  'tbp_pAT_choice_prob' : tbp_pAT_choice_probD,
                  'tbp_pAT_seq' : tbp_pAT_seqD,
                  'tf_sites' : tf_siteD,
                  'pAT_sites' : polyAT_siteD,
                  'tf_loci' : tf_lociD,
                  'pAT_loci' : polyAT_lociD,
                  'tf_add' : tf_addD,
                  'pAT_add' : pAT_addD,
                  'pAT_seqs' : pAT_seqD,
                  'tf_choice_prob' : tf_choice_probD,
                  'pAT_choice_prob' : pAT_choice_probD,
                  'reb1' : reb1,
                  'reb1_site_choice' : reb1_site_choiceD,
                  'rap1' : rap1,
                  'rap1_site_choice' : rap1_site_choiceD,
                  'gcr1' : gcr1,
                  'gcr1_site_choice' : gcr1_site_choiceD,
                  'abf1' : abf1,
                  'abf1_site_choice' : abf1_site_choiceL,
                  'mcm1' : mcm1,
                  'mcm1_site_choice' : mcm1_site_choiceL,
                  'rsc3' : rsc3,
                  're' : reD,
                  'atg_erase' : atgD,
                  'nn_erase' : nnD,
                  'scars' : scarD,
                  'flanks' : flankD}

    return parameterD

def cell(workbook, worksheet, cellname) :

    book = open_workbook(workbook)
        
    sheet = book.sheet_by_name(worksheet)

    row = int(re.sub("\D", "", cellname))-1
    col = re.sub("\d", "", cellname)
    
    col_list = generate_colnames()

    ncol = col_list.index(col)

    value = sheet.cell(row, ncol).value

    return value

def generate_colnames():
    
    col1 = product('ABCDEFGHIJKLMNOPQRSTUVWXYZ',repeat=1)
    col2 = product('ABCDEFGHIJKLMNOPQRSTUVWXYZ',repeat=2)
    cs = []

    for x in col1 :
        c = ''.join(x)
        cs.append(c)
        
    for y in col2 :
        c = ''.join(y)
        cs.append(c)

    return cs

if __name__ == "__maine__" :
    maine()
