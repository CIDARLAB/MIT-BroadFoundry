def maine() :
    rev_comp(raw_input("Seq:"))

def rev_comp(seq) :

    compD = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    
    seq_list = list(seq)

    rev_seq_list = []
    count = 0
    
    for x, nuc in enumerate(seq_list) :

        count = count + 1
    
        comp = compD[nuc]
        seq_list[x] = comp

        rev_seq_list.append('')

    rev_i = count-1
    
    for x, nuc in enumerate(rev_seq_list) :

        rev = seq_list[rev_i-x]  
        rev_seq_list[x] = rev

    rev_comp = ''.join(rev_seq_list)

    return rev_comp

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
