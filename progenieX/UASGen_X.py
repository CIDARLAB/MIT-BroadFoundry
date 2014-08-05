from xlrd import *
from common_functions import *
from excel_functions import cell
from reverser import rev_comp
from SeqGen_5 import seqgen
from pro_analysis import promoter_analysis
from site_eraser import re_eraser
from uas_subs import elmsub
    
def maine() :
    uasgen(int(raw_input("Number (must be a multiple of 4):")))
    
def uasgen(uasnum):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'

    # Define names of strength levels and iterable lists
    strengths = ['VH', 'H', 'M', 'L']
    ATCG = ['A', 'T', 'C', 'G']
    uas_num = [1,2]
    
    # These are the Type IIs sites to add to the end of the sequence for cloning
    bbs1_F = {1 : cell(iB,gS,'G7'), 2 : cell(iB,gS,'G9')}
    bbs1_R = {1 : rev_comp(cell(iB,gS,'G8')), 2 : rev_comp(cell(iB,gS,'G10'))}
    
    # These are scars I will be using for all generated sequences
    # For UAS1, the left scar will be odd J and the right even J for each strength
    # For UAS2, the left scar will be A and the right odd J for each strength
    # These scars are written in the Excel file
    A = cell(iB,gS,'B5')
    B = cell(iB,gS, 'B6')
    J1 = cell(iB,gS,'B7')
    J2 = cell(iB,gS,'B8')
    J3 = cell(iB,gS,'B9')
    J4 = cell(iB,gS,'B10')

    uas_scar = {1 : [J2, J3], 2 : [A, J4]}

    # Since the generator will create VH, H, L and M for each uasnum, to output the
    # number of desired sequences, the input has to be divided by 4.
    uass = uasnum/4

    for y in uas_num :

        uasD = {'VH' : [], 'H' : [], 'M' : [], 'L' : []}
        count = 1

        # Clear motif substitution file
        clear('uas%(y)s_sub_record.txt' % {'y': y})

        # Begin sequence generation for each strength level
        for x in strengths :
            
            # Generate random sequences with nucleotide percentages common in promoters
            seqgen(int(cell(iB,uS,'A3')), int(cell(iB,uS,'B3')),
                   int(cell(iB,uS,'C3')), int(cell(iB,uS,'D3')),
                   float(cell(iB,uS,'E3')), int(cell(iB,uS,'F3')), uass)

            # Retrieves the generated sequences from the output file and
            # makes a list of the sequences without the names added in seqgen()
            sygenlist = lister('seqsgen.txt')
            # Eliminate the generic FASTA names applied by seqgen
            uaslist = [line for line in sygenlist if '>' not in line]
            
            # Run the elmsub() function for each sequence and eliminate TypeIIS sites that
            # may have appeared. Concatenate sequences with the appropriate scar sequence.
            for num, seq in enumerate(uaslist) :
                with open('uas%(y)s_sub_record.txt' % {'y': y}, 'a') as rec :
                    rec.write('%(c)s %(x)s\n' % {'c': count, 'x': x})
                count = count + 1
                uas_el = elmsub(seq, x, y, iB, uS)
                uas = re_eraser(uas_scar[y][0]+uas_el+uas_scar[y][1])
                uaslist[num] = bbs1_F[y]+uas+bbs1_R[y]

            # Save list to the UAS dictionary defined earlier, output in FASTA format,
            # and analyze the sequences with the analysis suite.
            uasD[x] = uaslist
            prefix = 'uas%(num)s' % {'num': y}
            fastaD_out(uasD, prefix)
            promoter_analysis(prefix)
            
if __name__ == "__maine__" :
    maine()
