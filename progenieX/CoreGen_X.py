from xlrd import *
from common_functions import *
from excel_functions import cell
from reverser import rev_comp
from SeqGen_5 import seqgen
from core_element_subber import *
from site_eraser import *
from pro_analysis import promoter_analysis

def maine() :
    coregen(int(raw_input("Number (must be a multiple of 4):")))

def coregen(cornum) :
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    # Define names of strength levels and iterable lists 
    strengths = ['VH', 'H', 'M', 'L']
    subpart = ['tbp', 'tss', 'utr']
    ATCG = ['A', 'T', 'C', 'G']
    

    syxD = {'VH': {'tbp': [], 'tss': [], 'utr': []},
            'H': {'tbp': [], 'tss': [], 'utr': []},
            'M': {'tbp': [], 'tss': [], 'utr': []},
            'L': {'tbp': [], 'tss': [], 'utr': []}}

    coreD = {'VH': [], 'H': [], 'M': [], 'L': []}

    # These are all of the scars I will be using for all generated sequences
    # For core promoters, the left scar will be even J and right will be B
    A = cell(iB,gS,'B5')
    B = cell(iB,gS, 'B6')
    J1 = cell(iB,gS,'B7')
    J2 = cell(iB,gS,'B8')
    J3 = cell(iB,gS,'B9')
    J4 = cell(iB,gS,'B10')
    
    core_scar = [J1, B]

    # These are the Type IIs sites to add to the end of the sequence for cloning
    bbs1_F = cell(iB,gS,'G5')
    bbs1_R = rev_comp(cell(iB,gS,'G6'))
    
    # Since the generator will create VH, H, L and M for each cornum, to output the
    # number of desired sequences, the input has to be divided by 4.
    cors = cornum/4

    # Associated fractions from Lubliner et al. 2013
    specD = {'VH': {'tbp': {'A': cell(iB,cS,'C3'), 'T': cell(iB,cS,'D3'),
                            'C': cell(iB,cS,'E3'), 'G': cell(iB,cS,'F3')},
                    'tss': {'A': cell(iB,cS,'C4'), 'T': cell(iB,cS,'D4'),
                            'C': cell(iB,cS,'E4'), 'G': cell(iB,cS,'F4')},
                    'utr': {'A': cell(iB,cS,'C5'), 'T': cell(iB,cS,'D5'),
                            'C': cell(iB,cS,'E5'), 'G': cell(iB,cS,'F5')}},
             'H': {'tbp': {'A': cell(iB,cS,'C6'), 'T': cell(iB,cS,'D6'),
                           'C': cell(iB,cS,'E6'), 'G': cell(iB,cS,'F6')},
                    'tss': {'A': cell(iB,cS,'C7'), 'T': cell(iB,cS,'D7'),
                            'C': cell(iB,cS,'E7'), 'G': cell(iB,cS,'F7')},
                    'utr': {'A': cell(iB,cS,'C8'), 'T': cell(iB,cS,'D8'),
                            'C': cell(iB,cS,'E8'), 'G': cell(iB,cS,'F8')}},
             'M': {'tbp': {'A': cell(iB,cS,'C9'), 'T': cell(iB,cS,'D9'),
                           'C': cell(iB,cS,'E9'), 'G': cell(iB,cS,'F9')},
                    'tss': {'A': cell(iB,cS,'C10'), 'T': cell(iB,cS,'D10'),
                            'C': cell(iB,cS,'E10'), 'G': cell(iB,cS,'F10')},
                    'utr': {'A': cell(iB,cS,'C11'), 'T': cell(iB,cS,'D11'),
                            'C': cell(iB,cS,'E11'), 'G': cell(iB,cS,'F11')}},
             'L': {'tbp': {'A': cell(iB,cS,'C12'), 'T': cell(iB,cS,'D12'),
                           'C': cell(iB,cS,'E12'), 'G': cell(iB,cS,'F12')},
                    'tss': {'A': cell(iB,cS,'C13'), 'T': cell(iB,cS,'D13'),
                            'C': cell(iB,cS,'E13'), 'G': cell(iB,cS,'F13')},
                    'utr': {'A': cell(iB,cS,'C14'), 'T': cell(iB,cS,'D14'),
                            'C': cell(iB,cS,'E14'), 'G': cell(iB,cS,'F14')}}}
    
    toleranceD = {'tbp': cell(iB,cS,'G3'),
                  'tss': cell(iB,cS,'G4'),
                  'utr': cell(iB,cS,'G5')}
    
    lengthD = {'tbp': int(cell(iB,cS,'H3')),
               'tss': int(cell(iB,cS,'H4')),
               'utr': int(cell(iB,cS,'H5'))}
    
    clear('core_sub_record.txt')
    
    for x in strengths :
        
        for y in subpart :
            
            with open('core_sub_record.txt', 'a') as rec:
                rec.write('\n%(x)s\n%(y)s\n' % {'x':x, 'y':y})

            # Generate 50 bp sequences for each region
            seqgen(specD[x][y]['A'], specD[x][y]['T'],
                   specD[x][y]['C'], specD[x][y]['G'],
                   toleranceD[y], lengthD[y], cors)

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
                count = 1
                if y is 'tbp':
                    for num, sy in enumerate(sylist) :
                        with open('core_sub_record.txt', 'a') as rec:
                            rec.write(' %(c)s\n' % {'c':count})
                        sy = tatasub(sy, iB, cS)
                        count = count + 1
                        sylist[num] = sy

            # This function writes the polyAT into the tbp region.
            # This arrangement allows for one more pdW, but without combinatorial
            # explosion in UAS1, and also allows nucleosome free region close to the
            # tss. 
            if y is 'tbp':
                count = 1
                for num, sy in enumerate(sylist) :
                    with open('core_sub_record.txt', 'a') as rec:
                        rec.write(' %(c)s\n' % {'c':count})
                        sy = polyAT_tbp_sub(sy, x, iB, cS)
                        
            # Since the consensus Kozak is also unlikely to appear at the end of the UTR,
            # I also wrote a function that substitutes it in.
            if y is 'utr':
                count = 1
                for num, sy in enumerate(sylist) :
                    with open('core_sub_record.txt', 'a') as rec:
                        rec.write(' %(c)s\n' % {'c':count})
                    sy = kozaksub(sy, x, iB, cS)
                    sy = atg_eraser(sy, x)
                    sy = nab_nrd_eraser(sy, x)
                    count = count + 1
                    sylist[num] = sy

            # Since there is data available about elements around the TSS from Lubliner
            # I also wrote a function that puts in TSS elements right at the end of the
            # tss element
            if y is 'tss' :
                count = 1
                for num, sy in enumerate(sylist) :
                    with open('core_sub_record.txt', 'a') as rec:
                        rec.write(' %(c)s\n' % {'c':count})
                    sy = tsssub(sy, x, iB, cS)
                    sy = nab_nrd_eraser(sy, x)
                    count = count + 1
                    sylist[num] = sy

            # This populates the different segments into one dictionary
            syxD[x][y] = sylist

    # This stitches together the subsegments into a full core promoter and adds TypeIIS sites.
    for x in strengths :
        
        coreD[x] = syxD[x]['tbp']
        
        for num, seq in enumerate(syxD[x][y]):
            
            core = re_eraser(core_scar[0]+
                             syxD[x]['tbp'][num]+
                             syxD[x]['tss'][num]+
                             syxD[x]['utr'][num]+
                             core_scar[1])
            
            coreD[x][num] = bbs1_F+core+bbs1_R

    # Name output data files, format sequence dictionary into FASTA, and analyze for motifs.
    prefix = 'core'   
    fastaD_out(coreD, prefix)    
    promoter_analysis(prefix)


if __name__ == "__maine__" :
    maine()
