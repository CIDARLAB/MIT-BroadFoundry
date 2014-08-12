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
    
    # Slice Locations - called "Insertion Locations" in ProGenie_Parameters
    core_sliceD = {'tata'  : cell(iB,cS,'C18'),
                   'tss'   : cell(iB,cS,'C21'),
                   'kozak' : cell(iB,cS,'C22')}
    
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
    parameterD = {'nuc_pct': nuc_pctD,
                  'tolerance' : toleranceD,
                  'length' : lengthD,
                  'core_slice' : core_sliceD,
                  'tf_sites' : tf_siteD,
                  'pAT_sites' : polyAT_siteD,
                  'tf_loci' : tf_lociD,
                  'pAT_loci' : polyAT_lociD,
                  'tf_add' : tf_addD,
                  'pAT_add' : pAT_addD,
                  'pAT_seqs' : pAT_seqD,
                  'tf_choice_prob' : tf_choice_probD,
                  'pAT_choice_prob' : pAT_choice_probD,
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
