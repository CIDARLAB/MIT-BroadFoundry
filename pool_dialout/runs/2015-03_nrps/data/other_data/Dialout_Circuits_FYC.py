'''
Created on Sep 16, 2013

@author: dbg
'''

import sys
import argparse

#Ben-specific library
from DNA import FastaLoad
from DNA import align

def main():
    parser = argparse.ArgumentParser(description="Generate reference data for cistron pools")
    
    args = parser.add_argument_group("Arguments")
    
    args.add_argument('--partsfile',      '-p', action="store", default=None, dest="PartsFile",     help='CSV file containing pool part composition')
    args.add_argument('--partSeqs',       '-f', action="store", default=None, dest="PartSeqFile",   help='FAStA file containing all parts')
    args.add_argument('--invprimers',     '-i', action="store", default=None, dest="InvPrimers",    help='"ACAGT,TCACHA" primers for inverse PCR')
    args.add_argument('--linker',         '-L', action="store", default=None, dest="Linker",        help='linker sequence')
    args.add_argument('--regex',          '-R', action="store_true", default=False,dest="RegEx",         help='output regular expressions')
    options = parser.parse_args()
    
    partsfile  = options.PartsFile
    seqfile    = options.PartSeqFile
    invprimers = options.InvPrimers.upper().split(',')
    linker     = options.Linker
    regexflag  = options.RegEx
   
    if not partsfile or not seqfile:
        print "USAGE: Dialout_Cistrons.py -p partsfile.csv -f seqfile.fa"
        exit(1)
    
    id = partsfile.split('/')[-1].replace('.csv','')

    seqD = FastaLoad(seqfile)
    partslists = []
    cdsSet = set()
    
    for line in open(partsfile).readlines():
        parts = [x for x in line.strip().split(',') if x != '']
        partslists.append(parts)
    
    list_of_partslists = compute_list_of_partslists(partslists)    
    
    #seqnamepairs will be replaced at the end of each iteration with a new list
    #that contains all combinations so far
    
    ans = []
    for partslists in list_of_partslists:
        seqnamepairs = [([],['cistron'])]
    
        for partslist in partslists:
            newseqnamepairs = []
            for seqlist,namelist in seqnamepairs:
                for part in partslist:
                    newseqlist = seqlist + [seqD[part]]
                    newnamelist = namelist + [part]
                    newseqnamepairs.append((newseqlist, newnamelist))
            seqnamepairs = newseqnamepairs
        
        ans.extend(seqnamepairs)
                        
    if '_2NOR' in id:
        tokidxs = [3, 5, 7, 10, 18]
    elif '_NOT' in id:
        tokidxs = [3, 5 , 8, 16]
    elif '_FYC' in id:
        tokidxs = [5, 8]
    
    if not invprimers:
        j = 0
        for seqlist,namelist in ans:
            #seqid = name.replace('cistron','%s_%d'%(id,i+1))
            toks = namelist
            seqid = '_'.join(toks[i] for i in tokidxs)
            seq = ''.join(seqlist)
            print '>%s\n%s'%(seqid,seq)
            j += 1
    else:
        j = 0
        for seqlist,namelist in ans:
            toks = namelist
            seqid = '_'.join(toks[i] for i in tokidxs)
            preseq = ''.join(seqlist).upper()
            
            try:
                fwds = [m[0] for m in align(invprimers[0].upper(),preseq) if m[1] == 'f']
                fwdpos = fwds[0]
            except:
                print j, invprimers[0]
                print preseq
                fwds = [m[0] for m in align(invprimers[0].upper(),preseq) if m[1] == 'f']
                print fwds
                fwdpos = fwds[0]
                
            revs = [m[0] for m in align(invprimers[1].upper(),preseq) if m[1] == 'r']
            revpos = revs[0]
            
            seq = preseq[fwdpos:] + linker + preseq[:revpos]
            if regexflag: 
                seq = string2regex(seq)
            print '>%s\n%s'%(seqid,seq)
            j += 1
            
import re     
def string2regex(s):
    laststart = 0
    s = s.upper()
    Nres = list(re.finditer('N+',s))
    newstring = ''
    for Nre in Nres:
        newstring = newstring + s[laststart:Nre.start()]
        Nlen = Nre.end() - Nre.start()
        newstring = newstring + '([ACGT]{%d})'%(Nlen)
        laststart = Nre.end()
    newstring = newstring + s[laststart:]
    return newstring
    
def test_string2regex():    
    s = 'actNNNNTTTNNNTATATACTACNNNNNNNNNNNNNNC'
    print s
    print string2regex(s)
    s = 'actNNNNTTTNNNTATATACTACNNNNNNNNNNNNNN'
    print s
    print string2regex(s)
    s = 'NNNNTTTNNNTATATACTACNNNNNNNNNNNNNNC'
    print s
    print string2regex(s)

def compute_list_of_partslists(partslists):
    #Generate new partslists, only allowing f w f and r w r
    #There are allowed:

    #   RegTC_5 Bbs1 1 <more stuff> 9 Bbs1R
    #   RegTC_5 Bbs1 r2 <more stuff> r1 Bbs1R
    
    #and these are disallowed:
    
    #   RegTC_5 Bbs1 1 <more stuff> r1 Bbs1R
    #   RegTC_5 Bbs1 r2 <more stuff> 2 Bbs1R
    lopl = []
    
    #only Forward:
    ans = []
    for partslist in partslists:
        if all([part.replace('r','').isdigit() for part in partslist]):  #Only matches scars
            ans.append([part for part in partslist if part.find('r')<0]) #Keep only forw
        else:
            ans.append(partslist[:])
    lopl.append(ans)
    
    #only reverse:
    ans = []
    for partslist in partslists:
        if all([part.replace('r','').isdigit() for part in partslist]):  #Only matches scars
            ans.append([part for part in partslist if part.find('r')>=0]) #Keep only forw
        else:
            ans.append(partslist[:])
    lopl.append(ans)
    
    #print count_combinations(lopl[0]), lopl[0]
    #print count_combinations(lopl[1]), lopl[1]
    return lopl

def count_combinations(lol):
    ans = 1
    for sublist in lol:
        ans = ans * len(sublist)
    return ans
            

if __name__ == '__main__': main()
