

def maine():
    flow_reader(raw_input('Output filename:'))

def flow_reader(filename):

    with open('input.txt', 'r') as f :
        
        all_text = f.read()

    chunks = all_text.split('	0		')

    samples = all_text.split('	0		')

    for n, chunk in enumerate(chunks) :
        samples[n] = chunk.split('\n')

    
    dataD = {}
    well_list = []
    for n, sample in enumerate(samples) :
        sampleD = {}
        gfp = []
        rfp = []
        bfp = []
        name = sample[0]
        #well = name[name.find('-')+6:name.find('.')]
        well_list.append(n)
        for x, line in enumerate(sample) :
            
            if 'GFP-A' in line :
                gfp.append(line[line.find('Mean')+5:len(line)-1])
                sampleD['GFP'] = gfp
                    
            if 'mRFP-A' in line :
                rfp.append(line[line.find('Mean')+5:len(line)-1])
                sampleD['RFP'] = rfp
                
            dataD[n] = sampleD
            
    clear(filename)
    print dataD
    with open(filename, 'a') as out:
        out.write('Sample,GFP,RFP,\n')

        for x, well in enumerate(well_list) :
            print dataD[well_list[x]]
            if len(dataD[well_list[x]]) is not 0 :
            
                if len(dataD[well_list[x]]['GFP']) < 2:
                    out.write('%(w)s,%(g)s,%(r)s\n'% {'w': well,
                                                      'g':dataD[well]['GFP'][0],
                                                      'r':dataD[well]['RFP'][0]})

    clear('input.txt')

def clear(filename) :

    f = open(filename, 'w')
    f.write('')
    f.close()
    f = open(filename, 'a')

    return f  

if __name__ == "__maine__" :
    maine()
