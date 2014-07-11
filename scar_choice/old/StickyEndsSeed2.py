
from DNA import revcomplement
from itertools import product
from random import shuffle

seed = [("GGAG", "CTCC"),
        ("AGGT", "ACCT"),
        ("GCTT", "AAGC"),
        ("CGCT", "AGCG"),
        ("TACT", "AGTA"),
        ("TGCC", "GGCA"),
        ("GCAA", "TTGC"),
        ("CAGA", "TCTG"),
        ("CCAT", "ATGG"),
        ("AATG", "CATT")]

def main():
    se_pairs = generate_sticky_end_pairs(4)

    print "Found", len(se_pairs),"unique pairs of sticky ends"

    #I'm too lazy to write an exhaustive search algorithm, 
    #so let's instead randomize and do a greedy search,
    #keeping track of the best set found so far ("maxset").
    maxlen = 0
    maxset = None
    for i in range(5000):
        paircopy = se_pairs[:]
        shuffle(paircopy)
    #comment out next line for fully random search without seed
    paircopy = seed + paircopy
    orset = find_orthogonal(paircopy,min_hamming=2)
    if len(orset) > maxlen:
        print "Improved size of orthogonal set from",maxlen,"to",len(orset)
        maxlen = len(orset)
        maxset = orset.copy()
            
    for t in sorted(list(maxset)): print t
    #for t in list(maxset): print t
    print "Size of largest orthogonal set found:",len(maxset)
    
def find_orthogonal(se_pairs,min_hamming=2):
    'Greedy search, starting with first pair of sticky ends from those available'
    sets = []
    s = set()
    s.add(se_pairs[0])
    for f,r in se_pairs:
        keep = True
    #check for and eliminate palindromic scars
    if f == r:
        keep = False
        #continue
    #check if current scar contains at least 1 G or C
    if "G" not in f and "C" not in f:
        keep = False
        #continue
    #check if current scar contains no more than 3 G's or C's
    if "A" not in f and "T" not in f:
        keep = False
        #continue
        for curr_f,curr_r in s:
            #print 'c,f,r, h(f), h(r) =', curr_f, f, r, hamming(curr_f,f), hamming(curr_f,r)
            if hamming(curr_f,f) < min_hamming or hamming(curr_f,r) < min_hamming:
                keep = False
        if keep: 
            s.add((f,r))
    return s

def hamming(s1,s2):
    'Compute the hamming distance between two strings'
    h = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            h += 1
    return h
 
    
def generate_sticky_end_pairs(n=4):
    '''
    Generate an exhaustive list of all sticky ends, paired with
    their reverse complements (e.g. "ACTT", "AAGT"), and then filter
    for unique pairs.  For length n=4, there are 136 pairs.
    '''
    forw = product('ACGT',repeat=n)
    tups = []
    singles = set()
    for t in forw:
        f = ''.join(t)
        r = revcomplement(f)
        if (f not in singles) and (r not in singles):            
            tups.append((f,r))
            singles.add(f)
            singles.add(r)
    return tups
    

if __name__ == "__main__":
    '''This kicks off the program'''
    main()
