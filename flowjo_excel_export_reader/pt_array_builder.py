

def maine():
    pt_array_build(raw_input('File:'), raw_input('Print set:'), raw_input('Type (pro, ter, combo):'))

def pt_array_build(filename, printset, qtype) :

    with open(filename) as f :
        inlist = f.readlines()

    prolist = [line[:line.find('\t')] for line in inlist]

    terlist = [line[line.find('\t')+1:line.find('\t',3)] for line in inlist]

    gfplist = [line[line.find('\t',3)+1:len(line)-1] for line in inlist]
    
    qtypeD = {'pro': 1, 'ter': 2, 'combo' : 3}
    
    print printset

    if qtypeD[qtype] is 1 :

        for n, pro in enumerate(prolist):
            
            if int(pro) is int(printset) :
                print pro, terlist[n], gfplist[n]
                
        
    if qtypeD[qtype] is 2 :
        for n, ter in enumerate(terlist):
            if int(ter) is int(printset) :
                print prolist[n], ter, gfplist[n]
            

    if qtypeD[qtype] is 3:
       combo = printset.split(',')
       for n, pro in enumerate(prolist):

            if int(pro) is int(combo[0]) :
                if int(terlist[n]) is int(combo[1]):
                    print gfplist[n]
    

if __name__ == "__maine__" :
    maine()
