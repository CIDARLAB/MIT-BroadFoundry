

def maine():
    flow_reader(raw_input('Input filename:'))

def flow_reader(filename):

    with open(filename, 'r') as f :
        
        flist = f.readlines()
        
        for n, line in enumerate(flist) :
            
            flist[n] = line.replace('\n', '')

    glist = [line[31:] for line in flist if 'GFP-A' in line]

    for n, line in enumerate(glist) :
        glist[n] = line.replace('\t', '')

    rlist = [line[31:] for line in flist if 'mRFP-A' in line]

    for n, line in enumerate(rlist) :
        rlist[n] = line.replace('\t', '')

    nomlist = [line for line in flist if '.fcs' in line]
    dotlist = [line for line in flist if '.fcs' in line]
    hyplist = [line for line in flist if '.fcs' in line]
    
    for n, line in enumerate(dotlist) :
        dotlist[n] = line.find('.')

    for n, line in enumerate(hyplist) :
        hyplist[n] = line.find('-')
        
    for n, line in enumerate(nomlist) :
        nomlist[n] = line[hyplist[n]+6:dotlist[n]]

    for n, line in enumerate(nomlist) :
        print nomlist[n], glist[n]


if __name__ == "__maine__" :
    maine()
