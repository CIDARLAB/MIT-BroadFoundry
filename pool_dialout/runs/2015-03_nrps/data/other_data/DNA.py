'''
Created on Apr 23, 2013

@author: dbg
'''
import re
from   string import maketrans


##############
# Utilities
##############

def FastaLoad(filename):
    D = {}
    f = open(filename)
    chunks = [c.strip() for c in f.read().split('>') if len(c) > 1]
    for chunk in chunks:
        lines = chunk.split('\n')
        name = lines[0].split()[0]
        seq  = ''.join([line.strip() for line in lines[1:]])
        D[name] = seq
        #print ">>%s<<%d>>\n[%s]"%(name,len(seq),seq)
    return D
    
def FastaPrint(fastaD,sortedkeys=[]):
    if not sortedkeys: sortedkeys = sorted(fastaD.keys())
    for k in sortedkeys:
        seq = fastaD[k]
        print ">%s %d"%(k,len(seq))
        i = 0
        while i < len(seq):
            print seq[i:i+70]
            i+= 70

revcomplement_memo = {'A':'T'}
revcompTBL = maketrans("AGCTagctWSKMYRnN", "TCGAtcgaWSMKTYnN")
def revcomplement(seq):
    """
    revcomplement(seq)
    A quick reverse-complement routine that memo-izes queries, understands
    IUPAC ambiguity codes, and preserves case.
    """
    global revcomplement_memo
    try:
        rc = revcomplement_memo[seq]
    except KeyError:
        _t = list(seq.translate(revcompTBL))
        _t.reverse()
        rc = ''.join(_t)
        revcomplement_memo[seq] = rc
        revcomplement_memo[rc]  = seq
    return(rc)

def reverse_complement(x): 
    '''Wrapper for Bio compatibility'''
    return revcomplement(x)

def align(primer,seq):
    rcprimer = revcomplement(primer)
    forw = [(m.start(),"f")             for m in re.finditer('(?=%s)'%primer, seq)]
    rev  = [(m.start()+len(primer),"r") for m in re.finditer('(?=%s)'%rcprimer, seq)]
    return sorted(forw+rev)
