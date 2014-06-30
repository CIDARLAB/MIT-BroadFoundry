from common_functions import *
from CoreGen_X import coregen
from UASGen_X import uasgen
from pro_analysis import promoter_analysis

def maine() :
    pro_genie(int(raw_input("Number of sequences (must be multiple of 4):")))

def pro_genie(num) :

    # Generate sequences using the coregen and uasgen functions
    coregen(num)
    uasgen(num)

    # Retrieve generated sequences from the output txt files.  Omit TypeIIS sites on ends
    # to enable stitching as they would be in a TypeIIS reaction.
    corelist = lister('coregen.txt')

    cores = [line[11:169] for line in corelist if '>' not in line]
    x = [line[9:10] for line in corelist if '>' in line]

    uas1list = lister('uas1gen.txt')

    uas1s = [line[11:165] for line in uas1list if '>' not in line]
    
    uas2list = lister('uas2gen.txt')

    uas2s = [line[11:165] for line in uas2list if '>' not in line]

    # Stitch together the list of sequences into promoters. Note that this is not random stitching,
    # rather sequences are combined in the order they are generated.
    pros = cores

    for n, pro in enumerate(pros) :
        pro = uas2s[n]+uas1s[n]+cores[n]
        pros[n] = pro

    # Name output data files, format sequence dictionary into FASTA, and analyze for motifs.
    prefix = 'pro'
    fasta_out(pros, prefix)
    promoter_analysis(prefix)
            

if __name__ == "__maine__" :
    maine()
