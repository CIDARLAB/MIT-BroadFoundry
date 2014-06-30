

def maine():
    utr_lister(raw_input('Promoter Name:'), raw_input('Terminator Name:'))

def utr_lister(pname, tname):

    with open('Pro-UTR-Seq.txt', 'r') as utr5 :

        utr5list = utr5.readlines()

    nom5list = [line[:line.find('>')] for line in utr5list]
    
    slice5list = [line[line.find('>')+1:line.find('-')] for line in utr5list]

    seq5list = [line[line.find('-')+1:] for line in utr5list]

    proD = {}

    for x, nom in enumerate(nom5list):
        if slice5list[x] is not '':
            seq5list[x] = seq5list[x][len(seq5list[x])-int(slice5list[x]):len(seq5list[x])-1]
            proD[nom5list[x]] = seq5list[x]

    
    with open('Ter-UTR-Seq.txt', 'r') as utr3 :

        utr3list = utr3.readlines()

    nom3list = [line[:line.find('>')] for line in utr3list]

    slice3list = [line[line.find('>')+1:line.find('-')] for line in utr3list]

    seq3list = [line[line.find('-')+1:] for line in utr3list]

    terD = {}

    for x, nom in enumerate(nom3list):
        if slice3list[x] is not '':
            seq3list[x] = seq3list[x][:int(slice3list[x])]
            terD[nom3list[x]] = seq3list[x]


    with open('vGFP.txt', 'r') as gfpseq :
        gfplist = gfpseq.readlines()
        gfp = gfplist[0].replace('\n', '')

    transcript = proD[pname]+gfp+terD[tname]
    print '''5'UTR: ''',proD[pname]
    print '''3'UTR(5->3): ''',terD[tname]
    print '''3'UTR(3->5): ''',rev(terD[tname])
    print transcript

def rev(seq) :

    seq_list = list(seq)

    rev_seq_list = []
    count = 0
    
    for x, nuc in enumerate(seq_list) :

        count = count + 1
        
        seq_list[x] = nuc

        rev_seq_list.append('')

    rev_i = count-1
    
    for x, nuc in enumerate(rev_seq_list) :

        rev = seq_list[rev_i-x]  
        rev_seq_list[x] = rev

    rev = ''.join(rev_seq_list)

    return rev
            
if __name__ == "__maine__" :
    maine()
